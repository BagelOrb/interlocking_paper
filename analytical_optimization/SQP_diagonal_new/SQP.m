clear; % otherwise changes to the script aren't loaded until you restart MATLAB
constants;
design_variables;
objective_sqp;
constraints_sqp;

syms wa wb L F
syms l1 l2 l3 l4
syms wa_k wb_k L_k F_k


x = [wa; wb; L; F];
lambda = [l1; l2; l3; l4];
obj_min = 999999;

opts = optimoptions(@quadprog,'Algorithm','trust-region-reflective');

n_iter = 100;
delta = 1000;
lambda_k = ones(4,1);

% Set active constraints
h1 = g3_2;
h2 = g4a;
h3 = g5a;
h4 = g5b;
h = [sym(h1); sym(h2); sym(h3); sym(h4)];

[A, W, dfdx] =  getAandWmatrix(f, h, x, lambda);
    
dfdx_k = subs(dfdx, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]);
h_k = subs(h, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]);
W_k = subs(W, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]);
A_k = subs(A, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]); 

x_k = [0.4; 0.9; 3.2; 2.5];

wa_k = x_k(1);
wb_k = x_k(2);
L_k = x_k(3);
F_k = x_k(4);

for n = 1:n_iter
    
    lambda_k = double(lambda_k + 1/delta*eval(h_k));
    W_lambda = subs(W_k, [l1 l2 l3 l4], [lambda_k(1), lambda_k(2), lambda_k(3), lambda_k(4)]);

    % lb = [wa_min - x_eval(1) , wb_min -  x_eval(2), L_min - x_eval(3), 99999]
    % ub = [w_max - x_eval(1), w_max - x_eval(2), L_max, -99999]
    
    
    d = x - x_k;
    W_eval = double(eval(W_lambda));
    dfdx_eval = double(eval(dfdx_k));
    A_eval = double(eval(A_k));
    h_eval = double(eval(h_k));
    
%     f_lin = @(wa, wb, L, F) subs(f) + dot(subs(dfdx), d) + 1/2*d.' * subs(W) * d;
%     h1_lin = @(wa, wb, L, F) subs(h1) + dot(A_eval(:,1), d);
%     h2_lin = @(wa, wb, L, F) subs(h2) + dot(A_eval(:,2), d);
%     h3_lin = @(wa, wb, L, F) subs(h3) + dot(A_eval(:,3), d);
%     h4_lin = @(wa, wb, L, F) subs(h4) + dot(A_eval(:,4), d);
%     
%     h_lin = [sym(h1_lin); sym(h2_lin); sym(h3_lin); sym(h4_lin)];
%     
%     [A_lin, W_lin, dfdx_lin] =  getAandWmatrix(f_lin, h_lin, x, lambda);
%     
%     W_lin_eval = double(subs(W_lin));
%     dfdx_lin_eval = double(subs(dfdx_lin));
%     A_lin_eval = double(subs(A_lin));
%     h_lin_eval = double(subs(h_lin));
    
    dx = quadprog(W_eval, dfdx_eval.', A_eval, -h_eval);
 
    x_k = x_k + dx
    
    ev = eig(W_eval)
    isposdef = all(ev>-0);
    
    wa_k = x_k(1);
    wb_k = x_k(2);
    L_k = x_k(3);
    F_k = x_k(4);
    
    obj  = f(wa_k, wb_k, L_k, F_k);
    
    
     if isposdef
        obj_min  = obj;
        break
    end
end


fprintf('The minimum objective of %f with a max nominal stress of %f is reached for: \n', obj_min, 1 / obj_min)
fprintf('w_a* = %f \n', wa);
fprintf('w_b* = %f \n', wb);
fprintf('L* = %f \n', L);
fprintf('F* = %f \n', F);


fprintf('And occurs for the following constraint values: \n');
fprintf('g1a = %f \n', g1a(wa, wb, L, F));
fprintf('g1b = %f \n', g1b(wa, wb, L, F));
fprintf('g2 = %f \n', g2(wa, wb, L, F));
fprintf('g3_1 = %f \n', g3_1(wa, wb, L, F));
fprintf('g3_2 = %f \n', g3_2(wa, wb, L, F));
fprintf('g4a = %f \n', g4a(wa, wb, L, F));
fprintf('g4b = %f \n', g4b(wa, wb, L, F));
fprintf('g5a = %f \n', g5a(wa, wb, L, F));
fprintf('g5b = %f \n', g5b(wa, wb, L, F));



