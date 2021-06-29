clear all; % otherwise changes to the script aren't loaded until you restart MATLAB
constants;
design_variables;
objective_sqp;
constraints_sqp;

% TODO:
% DONE - merge x_k x_k x_k and x_0 into x_k 
% DONE - make x_k into a list, remove all x_k(1), x_k(2), etc
% DONE - fix stopping criteria
% - define 2d subproblems
% - make 2d plots of f and include constraints
% - implement active set strategy
%     > use submatrices of A
%     > compute W based on full list of matrices of 2nd
%       order derivatives of ddhdxx

syms wa wb L F
syms l1 l2 l3 l4 l5 l6 l7 l8 l9 l10

% Define starting point
x = [wa; wb; L; F];
lambdas = [l1; l2; l3; l4; l5; l6; l7; l8; l9; l10]; % is needed for symbolic different

% Set index of active constraints
h_idx = [];
lambda_k = [];

% Set number of iterations
Niter = 25;

f_sym = sym(f);

ng = length(gs);
nh = length(h_idx);
nx = length(x);

nhistory = 6;
obj_history = [99999:99999+nhistory].';
h_idx_history = cell(nhistory,1);
cycling_break = 0;

x_k = ones(1,nx);
lambdas_k = ones(1,ng);

delta = 100;
options = optimset('Display', 'off');

% Obtain second order Taylor series approximation
disp('Approximating Taylor series...');
g = cellfun(@(g_) sym(g_), gs);
h = g(h_idx);
h_lin = eval(h);


g_eval = ones(ng,1);

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

for p = 1:Niter
    % Determine active constraints
    g_k = eval(subs(g, x.', x_k));
    lambdas_k(h_idx) = lambda_k;
    
    % cycle detection and resolution
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
        h_idx = [indices(g_k > .000001)];
    end
    nh = length(h_idx);
    h = g(h_idx);
    lambda = lambdas(h_idx); % change subset of lambdas to match subset of active constraints
    lambda_k = lambdas_k(h_idx);
    
    % Calculate Lagrangian multipliers
    h_k = g_k(h_idx);
    lambda_k = double(lambda_k + 1/delta*double(h_k));
    
    % Compute Hessian and Jacobian matrices for SQP
    ddhdxx = ddgdxx(:,:,h_idx);
    W = ddfdxx;
    for i = 1:nh
        W = W + lambda_k(i) * ddhdxx(:,:,i);
    end
    A = dgdx(:,h_idx).';
    
    % Evaluate matrices at the current point
    dfdx_k = subs(dfdx, x.', x_k);
    W_k = subs(W, x.', x_k);
    A_k = subs(A, x.', x_k); 
    
    % Evaluate matrices with lambda
    W_eval = double(subs(W_k, lambda.', lambda_k));
    dfdx_eval = double(subs(dfdx_k, lambda.', lambda_k));
    A_eval = double(A_k);
    h_eval = double(h_k);

    % Obtain update step
    [dx, sqp_obj, exitflag, output, lambda_next] = quadprog(W_eval, dfdx_eval.', [], [], A_eval, -h_eval, [], [], [], options);
    x_k = x_k + dx.';
    
    if all(abs(dx) < 10^(-10))
        break;
    end

    % Evaluate eigenvalues
    ev = eig(W_eval);
    isposdef = all(ev> -10^(-14));

    % Compute objective
    obj  = eval(subs(f_sym, x.', x_k));
    g_eval = eval(subs(g, x.', x_k));
    fprintf("%i: objective: %.5f,\t #constraints: %i,\t highest constraint: %.3f\n", p, obj, length(h_idx), max(g_eval));
    
    % Record history
    obj_history = [obj_history(2:nhistory); obj];
    h_idx_history(1:nhistory-1,1) = h_idx_history(2:nhistory,1);
    h_idx_history(nhistory,1) = { h_idx };
end

fprintf('\n The minimum objective of %f with a max nominal stress of %f is reached for: \n', obj, 1 / obj)
for i = 1:nx
    fprintf('%s = %f \n', string(x(i)), x_k(i));
end

fprintf('And occurs for the following constraint values: \n');
for i = 1:ng
    fprintf('%s: %.3f \n', g_names(i), g_eval(i));
end