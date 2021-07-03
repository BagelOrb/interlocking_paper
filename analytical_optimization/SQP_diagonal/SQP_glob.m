clear all; % otherwise changes to the script aren't loaded until you restart MATLAB

straight_case_2var;
% diagonal_case_2var;
% straight_case;
% diagonal_case;


% Set optimization parameters
Niter = 1000;

show_every_iteration = false;

get_plot = 1;       % Turn on to obtain plot
syms l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 l12 l13

% Define starting point
lambdas = [l1; l2; l3; l4; l5; l6; l7; l8; l9; l10; l11; l12; l13]; % is needed for symbolic different

options = optimset('Display', 'off');
thres = 10^(-10);

% history for cycling detection
nhistory = 6;
obj_history = [99999:99999+nhistory].';
h_idx_history = cell(nhistory,1);
h_max_history = cell(nhistory,1);
x_history = [x_k];
cycling_break = 0;

% general initialization
h_idx = [];

ng = length(g);
nh = length(h_idx);
nx = length(x);

lambdas_k = ones(1,ng);
lambda_k = [];

best_x = ones(1, nx);
best_obj = 99999;

% Obtain second order derivatives
disp('Approximating derivatives...');
tic;

dfdx = [];
for i = 1:nx
    dfdx = [dfdx ; diff(f, x(i))];
end

ddfdxx = diff2(sym(f), x);
ddgdxx = sym(zeros(nx, nx, ng));
dgdx = sym(zeros(nx, ng));
for i = 1:ng
    [ddgdxx_i, dgdx_i] = diff2(g(i), x);
    ddgdxx(:,:,i) = ddgdxx_i;
    dgdx(:,i) = dgdx_i;
end

g_f = matlabFunction(g, 'Vars', x);
f_f = matlabFunction(f, 'Vars', x);
dfdx_f = matlabFunction(dfdx, 'Vars', x);
ddfdxx_f = matlabFunction(ddfdxx, 'Vars', x);
ddgdxx_f = matlabFunction(ddgdxx, 'Vars', x);
dgdx_f = matlabFunction(dgdx, 'Vars', x);

toc; % print time elapsed for differentiation

tic;
disp('Start iteration...');

for p = 1:Niter

    x_kc = num2cell(x_k, 1);
    dfdx_k = dfdx_f(x_kc{:});
    dgdx_k = dgdx_f(x_kc{:});
    dgdx_k(isnan(dgdx_k)) = 0;
    dgdx_k(isinf(dgdx_k)) = 0;

    % Determine active constraints
    g_k = g_f(x_kc{:});
    lambdas_k(h_idx) = lambda_k; % save lambdas associated with old h_idx


    % Check KKT conditions
    if ~ any(isnan(dgdx_k))
        if p > 1
            mu = lambdas_k;
            constraints_satisfied = all(g_k < thres);
            positive_mu = all(lambdas_k > -thres);
            inactive_or_satisfied = all(abs(lambdas_k .* g_k) < thres);

            if constraints_satisfied ...
                    && positive_mu...
                    && inactive_or_satisfied
                fprintf("Constraints satisfied.\n");
                break;
            end
        end
    end

    % Cycle detection and resolution
    if p > nhistory && ~ cycling_break ...
        && ~ isequal(cell2mat(h_idx_history(1)), cell2mat(h_idx_history(2)))...
        && isequal(cell2mat(h_idx_history(1)), cell2mat(h_idx_history(3)))...
        && isequal(cell2mat(h_idx_history(3)), cell2mat(h_idx_history(5)))...
        && isequal(cell2mat(h_idx_history(2)), cell2mat(h_idx_history(4)))...
        && isequal(cell2mat(h_idx_history(4)), cell2mat(h_idx_history(6)))
        fprintf("Cycling detected! Take a break and choose the one with more active constraints...\n");
        cycling_break = 1;
        h_idx = cell2mat(h_idx_history(1));
        if length(h_idx) < length(cell2mat(h_idx_history(2)))
            h_idx = cell2mat(h_idx_history(2));
        end
    end

    if ~ cycling_break
        indices = [1:ng]; % set all violated constraints as active 
        mu = lambdas_k;
        h_idx_before = h_idx;
        violateds = [indices(g_k > thres)];
        h_idx = violateds;

        if ~ isequal(h_idx_before, h_idx)
            lambdas_k = ones(1, ng);
        end

    else
        cycling_break = cycling_break + 1;
    end
    if cycling_break > 6
        cycling_break = 0;
    end

    nh = length(h_idx);
    h = g(h_idx);
    lambda = lambdas(h_idx); % change subset of lambdas to match subset of active constraints
    lambda_k = lambdas_k(h_idx);

    % Calculate Lagrangian multipliers
    h_k = g_k(h_idx);

    % Compute Hessian and Jacobian matrices for SQP
    ddgdxx_k = ddgdxx_f(x_kc{:});
    ddhdxx_k = ddgdxx_k(:,:,h_idx);
    W_k = ddfdxx_f(x_kc{:});
    for i = 1:nh
        W_k = W_k + lambda_k(i) * ddhdxx_k(:,:,i);
    end
    A_k = dgdx_k(:,h_idx).';

    % Obtain update step
    [dx, sqp_obj, exitflag, output, lambda_next] = quadprog(W_k, dfdx_k.', [], [], A_k, -h_k, [], [], [], options);

    if exitflag == -2
        fprintf("Problem non-convex! Stopping execution!\n");

        % Go one step back
        x_k = x_history(p - 1,:);
        x_history = x_history(1:end-1,:);
        obj  = f_f(x_kc{:});
        g_k = g_f(x_kc{:});
        break;
    else

        [lambda_new, rank] = linsolve(A_k.', - W_k * dx - dfdx_k);
        lambda_k(~ isinf(lambda_new)) = lambda_new(~ isinf(lambda_new));
    end

    employed_move_limits = false;
    lll = sum(dx .* dx);
    if lll > move_limit^2
        employed_move_limits = true;
        dx = dx * move_limit / sqrt(lll);
    end

    x_k = x_k + dx.';
    x_history = [x_history; x_k];   

    % Evaluate eigenvalues
    ev = eig(W_k);
    isposdef = all(ev> -thres);

    % Compute objective
    x_kc = num2cell(x_k, 1);
    obj  = f_f(x_kc{:});
    g_k = g_f(x_kc{:});

    if obj < best_obj && all(g_k < thres)
        best_obj = obj;
        best_x = x_k;
        best_lambda = lambda_k;
    end

    if show_every_iteration
        fprintf("%i: objective: %.5f,\t constraints: %s,\t highest constraint: %.3f,\t move limits: %i\n", p, obj, num2str(h_idx), max(g_k), employed_move_limits);
    end

    if max(g_k) < thres && max(g_k) > -0.01 && any(round(obj_history(:), 5) == round(obj,5)) && ismember(round(x_k,5), round(x_history,5), 'rows') && length(h_idx) == length(x)
        fprintf("Cycling inside the same loop, stop!\n")
        break;
    end

    % Record history
    obj_history = [obj_history(2:nhistory); obj];
    h_idx_history(1:nhistory-1,1) = h_idx_history(2:nhistory,1);
    h_idx_history(nhistory,1) = { h_idx };

    h_max_history(1:nhistory-1,1) = h_max_history(2:nhistory,1);
    h_max_history(nhistory,1) = { max(g_k) };

    if all(abs(dx) < thres)
        fprintf("Not budging anymore! dx too small.\n")
        break;
    end

    g_k_other = g_k;
    g_k_other(h_idx) = 0;
%     if max(g_k_other) > thres && cycling_break 
%         fprintf("Another constraint violated, stop!\n")
% 
%         % Go one step back
%         x_k = x_k - dx.';
%         obj  = f_f(x_kc{:});
%         g_k = g_f(x_kc{:});
%         x_history = x_history(1:end-1,:);
%         break
%     end

    if p == Niter
        fprintf("Reaching max iterations.\n");
    end
end
x_kc = num2cell(x_k, 1);
obj  = f_f(x_kc{:});

toc; % print time elapsed for iterations

if obj > best_obj || any(g_k > thres)
    fprintf("Taking best x...\n");
    obj = best_obj;
    x_k = best_x;
    lambda_k = best_lambda;
    x_kc = num2cell(x_k, 1);
    g_k = g_f(x_kc{:});
    dfdx_k = dfdx_f(x_kc{:});
    dgdx_k = dgdx_f(x_kc{:});
    dgdx_k(isnan(dgdx_k)) = 0;
    dgdx_k(isinf(dgdx_k)) = 0;
end

% Check KKT conditions
if ~ any(isnan(dgdx_k))
    if p > 1
        [mu, r] = linsolve(dgdx_k, -dfdx_k);
        constraints_satisfied = all(g_k < thres);
        positive_mu = all(mu > -thres);
        inactive_or_satisfied = all(abs(mu.' .* g_k) < thres);
        if constraints_satisfied ...
                && positive_mu...
                && inactive_or_satisfied
            fprintf("Constraints satisfied.\n");
        else
            fprintf("Constrains violated!\n");
        end
    end
end

fprintf('\n The minimum objective of %f with a max nominal stress of %f is reached for: \n', obj, 1 / obj)
for i = 1:nx
    fprintf('%s = %f \n', string(x(i)), x_k(i));
end

fprintf('And occurs for the following constraint values: \n');
for i = 1:ng
    fprintf('%s: %.3f \n', g_names(i), g_k(i));
end


if get_plot && nx == 2
    fprintf('Plotting...');
    x1_range = max(x_history(:,1)) - min(x_history(:,1));
    x2_range = max(x_history(:,2)) - min(x_history(:,2));
    x1_array = linspace(max(0.1, min(x_history(:,1)) - 0.1*x1_range), .1*x1_range+max(x_history(:,1)), 50); 
    x2_array = linspace(max(0.5, min(x_history(:,2)) - 0.1*x2_range), .1*x2_range+max(x_history(:,2)), 50);
    contourplots(x, x_history, x1_array, x2_array, f_f, g_f, g_names)
end
fprintf('Done!');

% Logarithmic sensitivities
for i = 1:size(dgdx_k, 1)
    for j = 1:size(dgdx_k, 2)
        dgdx_L_ij = x(i)/ g(j) * dgdx(i,j);
        dgdx_L_ij_f = matlabFunction(dgdx_L_ij, 'Vars', x);
        dgdx_L(i,j) = round(dgdx_L_ij_f(x_kc{:}),3);
    end
end

for i = 1:size(dfdx_k, 1)
        dfdx_L(i) = round(x_k(i)/best_obj * dfdx_k(i),3);    
end
dfdx_L = dfdx_L.';