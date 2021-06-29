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
syms l1 l2 l3 l4

% Define starting point
x = [wa; wb; L; F];
lambda = [l1; l2; l3; l4]; % is needed for symbolic different

% Set index of active constraints
h_idx = [1; 5; 8; 9];

% Set number of iterations
Niter = 100;       % Max number of iterations

f_sym = sym(f);

ng = length(gs);
nh = length(h_idx);
nx = length(x);

x_k = ones(1,nx);
lambda_k = ones(1,nx);

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
    q = 1;
    
    for g_con = g
        g_eval(q) = eval(subs(g_con, x.', x_k));
        q = q + 1;
    end
    
    % If the optimum is found -> stop iterating
    [g_active, idx]  = maxk(g_eval,nx);
    
    % Compute Hessian and Jacobian matrices for SQP
    h = g(h_idx);
    ddhdxx = ddgdxx(:,:,h_idx);
    W = ddfdxx;
    for i = 1:nh
        W = W + lambda_k(i) * ddhdxx(:,:,i);
    end
    A = dgdx(:,h_idx).';
    
    % Evaluate matrices at the current point
    dfdx_k = subs(dfdx, x.', x_k);
    h_k = subs(h, x.', x_k);
    W_k = subs(W, x.', x_k);
    A_k = subs(A, x.', x_k); 
    
    % Calculate Lagrangian multipliers
    lambda_k = double(lambda_k + 1/delta*double(h_k));
    % any update for lambda seems to be good; something is wrong!
    
    % Evaluate matrices with lambda
    W_eval = double(subs(W_k, lambda.', lambda_k));
    dfdx_eval = double(subs(dfdx_k, lambda.', lambda_k));
    A_eval = double(A_k);
    h_eval = double(h_k);

    % Obtain update step
    [dx, sqp_obj, exitflag, output, lambda_next] = quadprog(W_eval, dfdx_eval.', [], [], A_eval, -h_eval, [], [], [], options);
    x_k = x_k + dx.';
    %lambda_k = -lambda_next.eqlin.';
    
    if all(abs(dx) < 10^(-10))
        break;
    end

    % Evaluate eigenvalues
    ev = eig(W_eval);
    isposdef = all(ev> -10^(-14));

    % Compute objective
    obj  = eval(subs(f_sym, x.', x_k));     
    fprintf("%i: objective: %.3f,\t highest constraint: %.3f\n", p, obj, max(g_eval));
end

fprintf('\n The minimum objective of %f with a max nominal stress of %f is reached for: \n', obj, 1 / obj)
for i = 1:nx
    fprintf('%s = %f \n', string(x(i)), x_k(i));
end

fprintf('And occurs for the following constraint values: \n');
for i = 1:ng
    fprintf('%s: %.3f \n', g_names(i), g_eval(i));
end