clear all; % otherwise changes to the script aren't loaded until you restart MATLAB
constants;
design_variables;
objective_sqp;
constraints_sqp;

% Define starting point
x_0 = [1; 1; 1; 1];

% Set index of active constraints
h_idx = [1; 5; 8; 9];

% Set number of iterations
Niter = 1000;       % Max number of iterations

syms wa wb L F
syms l1 l2 l3 l4
syms wa_k wb_k L_k F_k
syms wa_iter wb_iter L_iter F_iter

x = [wa; wb; L; F];
r = length(x_0);
lambda = [l1; l2; l3; l4];
lambda_k = zeros(r,1);
obj_min = 999999;

wa_best = x_0(1);
wb_best = x_0(2);
L_best = x_0(3);
F_best = x_0(4);

delta = 100;
options = optimset('Display', 'off');

% Obtain second order Taylor series approximation
disp('Approximating Taylor series...');
g = cellfun(@(g_) sym(g_), gs);
g_lin = cellfun(@(g_) taylor2(g_, wa_k, wb_k, L_k, F_k), gs);

f_quad = taylor2(f, wa_k, wb_k, L_k, F_k);

g_eval = ones(length(g),1);

g_lin_eval = ones(length(g_lin),1);


% h = [g(h_idx(1)); g(h_idx(2)); g(h_idx(3)); g(h_idx(4))];
% [A_real, W_real, dfdx_real] =  getAandWmatrix(sym(f), h, x, lambda);

for p = 1:Niter
    disp(['Start iteration ', num2str(p)]);
    
    if p == 1
        x_k = x_0;
        wa_k = x_k(1);
        wb_k = x_k(2);
        L_k = x_k(3);
        F_k = x_k(4);
          
    else 
        x_k = x_best;
        wa_k = wa_best;
        wb_k = wb_best;
        L_k = L_best;
        F_k = F_best;
        clear wa_iter wb_iter L_iter F_iter
        syms wa_iter wb_iter L_iter F_iter
        
    end
    
    % Determine active constraints
    q = 1;
    for g_con = g
                    g_eval(q) = eval(subs(g_con, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best]));
                    q = q + 1;
    end
    
    % If the optimum is found -> stop iterating
    [g_active, idx]  = maxk(g_eval,r);
    if all(g_active >= -10^(-14)) && all(g_eval <= 10^(-14))
        break;
    end
    
    h_lin = [g_lin(h_idx(1)); g_lin(h_idx(2)); g_lin(h_idx(3)); g_lin(h_idx(4))];
    h_lin = eval(h_lin);

    % Find Hessian and Jacobian matrices for SQP
    [A, W, dfdx] =  getAandWmatrix(f_quad, h_lin, x, lambda);

    dfdx_k = subs(dfdx, [wa, wb, L, F], [wa_iter, wb_iter, L_iter, F_iter]);
    h_k = subs(h_lin, [wa, wb, L, F], [wa_iter, wb_iter, L_iter, F_iter]);
    W_k = subs(W, [wa, wb, L, F], [wa_iter, wb_iter, L_iter, F_iter]);
    A_k = subs(A, [wa, wb, L, F], [wa_iter, wb_iter, L_iter, F_iter]); 

    x_iter = x_k;
    wa_iter = x_k(1);
    wb_iter = x_k(2);
    L_iter = x_k(3);
    F_iter = x_k(4);
    
    obj_min = 9999;
        
    % Calculate Lagrangian multipliers
    lambda_k = double(lambda_k + 1/delta*eval(h_k));
    W_lambda = subs(W_k, [l1 l2 l3 l4], [lambda_k(1), lambda_k(2), lambda_k(3), lambda_k(4)]);

    W_eval = double(eval(W_lambda));
    dfdx_eval = double(eval(subs(dfdx_k, [l1 l2 l3 l4], [lambda_k(1), lambda_k(2), lambda_k(3), lambda_k(4)])));
    A_eval = double(eval(A_k));
    h_eval = double(eval(h_k));

    % Obtain update step
    dx = quadprog(W_eval, dfdx_eval.', [], [], A_eval, -h_eval, [], [], [], options);
    x_iter = x_iter + dx;

    % Evaluate eigenvalues
    ev = eig(W_eval);
    isposdef = all(ev> -10^(-14));

    wa_iter = x_iter(1);
    wb_iter = x_iter(2);
    L_iter = x_iter(3);
    F_iter = x_iter(4);

    % Conpute objective
    obj  = f(wa_iter, wb_iter, L_iter, F_iter);

    if obj < obj_min
        obj_min  = obj;
        wa_best = wa_iter;
        wb_best = wb_iter;
        L_best = L_iter;
        F_best = F_iter;
        x_best = [wa_best; wb_best; L_best; F_best];   
        % lambda_opt = lambda_k;
         
        
        
    end
    
%      Hessian of full obj function without Taylor approximation
%     W_lambda_real = subs(W_real, [l1 l2 l3 l4], [lambda_opt(1), lambda_opt(2), lambda_opt(3), lambda_opt(4)]);
%     W_real_eval = double(subs(W_lambda_real, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best]));
%     ev_real = eig(W_real_eval);
            
    
end

fprintf('\n The minimum objective of %f with a max nominal stress of %f is reached for: \n', obj_min, 1 / obj_min)
fprintf('w_a* = %f \n', wa_best);
fprintf('w_b* = %f \n', wb_best);
fprintf('L* = %f \n', L_best);
fprintf('F* = %f \n', F_best);

fprintf('And occurs for the following constraint values: \n');
fprintf('g1a = %f \n', eval(subs(sym(g1a), [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g1b = %f \n', eval(subs(sym(g1b), [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g2 = %f \n', eval(subs(sym(g2), [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g3_1 = %f \n', eval(subs(sym(g3_1), [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g3_2 = %f \n', eval(subs(sym(g3_2), [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g4a = %f \n', eval(subs(sym(g4a), [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g4b = %f \n', eval(subs(sym(g4b), [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g5a = %f \n', eval(subs(sym(g5a), [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g5b = %f \n', eval(subs(sym(g5b), [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g6 = %f \n', eval(subs(sym(g6), [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));