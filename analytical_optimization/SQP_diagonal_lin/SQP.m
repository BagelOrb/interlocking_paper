clear; % otherwise changes to the script aren't loaded until you restart MATLAB
constants;
design_variables;
objective_sqp;
constraints_sqp;

syms wa wb L F
syms l1 l2 l3 l4
syms wa_k wb_k L_k F_k
syms wa_iter wb_iter L_iter F_iter

x = [wa; wb; L; F];
lambda = [l1; l2; l3; l4];
obj_min = 999999;

%opts = optimoptions(@quadprog,'Algorithm','trust-region-reflective');

n_iter = 100;
delta = 100;
lambda_k = ones(4,1);

% Set active constraints
h1 = taylor2(g3_2, wa_k, wb_k, L_k, F_k);
h2 = taylor2(g4a, wa_k, wb_k, L_k, F_k);
h3 = taylor2(g5a, wa_k, wb_k, L_k, F_k);
h4 = taylor2(g5b, wa_k, wb_k, L_k, F_k);
h = [sym(h1); sym(h2); sym(h3); sym(h4)];


g1a_lin =  taylor2(g1a, wa_k, wb_k, L_k, F_k);
g1b_lin = taylor2(g1b, wa_k, wb_k, L_k, F_k);
g2_lin = taylor2(g2, wa_k, wb_k, L_k, F_k);
g3_1_lin = taylor2(g3_1, wa_k, wb_k, L_k, F_k);
g3_2_lin = taylor2(g3_2, wa_k, wb_k, L_k, F_k);
g4a_lin = taylor2(g4a, wa_k, wb_k, L_k, F_k);
g4b_lin = taylor2(g4b, wa_k, wb_k, L_k, F_k);
g5a_lin = taylor2(g5a, wa_k, wb_k, L_k, F_k);
g5b_lin = taylor2(g5b, wa_k, wb_k, L_k, F_k);
f_lin = taylor2(f, wa_k, wb_k, L_k, F_k);


x_k = [0.6; 1.4; 3; 3.2222];

wa_k = x_k(1);
wb_k = x_k(2);
L_k = x_k(3);
F_k = x_k(4);

f_quad = taylor2(f, wa_k, wb_k, L_k, F_k);
h_lin = subs(h);

[A, W, dfdx] =  getAandWmatrix(f_quad, h_lin, x, lambda);


dfdx_k = subs(dfdx, [wa, wb, L, F], [wa_iter, wb_iter, L_iter, F_iter]);
h_k = subs(h_lin, [wa, wb, L, F], [wa_iter, wb_iter, L_iter, F_iter]);
W_k = subs(W, [wa, wb, L, F], [wa_iter, wb_iter, L_iter, F_iter]);
A_k = subs(A, [wa, wb, L, F], [wa_iter, wb_iter, L_iter, F_iter]); 

x_iter = [0.6; 1.4; 3; 3.2222];
wa_iter = x_iter(1);
wb_iter = x_iter(2);
L_iter = x_iter(3);
F_iter = x_iter(4);

for n = 1:n_iter
    
    lambda_k = double(lambda_k + 1/delta*eval(h_k));
    W_lambda = subs(W_k, [l1 l2 l3 l4], [lambda_k(1), lambda_k(2), lambda_k(3), lambda_k(4)]);
     
%   lb = [wa_min - x_iter(1) , wb_min -  x_iter(2), L_min - x_iter(3), 99999]
%   ub = [w_max - x_iter(1), w_max - x_iter(2), L_max, -99999]
  
    d = x - x_k;
    W_eval = double(eval(W_lambda));
    dfdx_eval = double(eval(dfdx_k));
    A_eval = double(eval(A_k));
    h_eval = double(eval(h_k));
    
    dx = quadprog(W_eval, dfdx_eval.', [], [], A_eval, -h_eval);
 
    x_iter = x_iter + dx;
    
    ev = eig(W_eval)
    isposdef = all(ev> -10^(-14));
    
    wa_iter = x_iter(1);
    wb_iter = x_iter(2);
    L_iter = x_iter(3);
    F_iter = x_iter(4);
    
    obj  = f(wa_iter, wb_iter, L_iter, F_iter);
    
    isposdef = 1;
     if isposdef
        obj_min  = obj;
        wa_best = wa_iter;
        wb_best = wb_iter;
        L_best = L_iter;
        F_best = F_iter;
    end
end


fprintf('The minimum objective of %f with a max nominal stress of %f is reached for: \n', obj_min, 1 / obj_min)
fprintf('w_a* = %f \n', wa_best);
fprintf('w_b* = %f \n', wb_best);
fprintf('L* = %f \n', L_best);
fprintf('F* = %f \n', F_best);


fprintf('And occurs for the following constraint values: \n');
fprintf('g1a = %f \n', eval(subs(g1a_lin, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g1b = %f \n', eval(subs(g1b_lin, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g2 = %f \n', eval(subs(g2_lin, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g3_1 = %f \n', eval(subs(g3_1_lin, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g3_2 = %f \n', eval(subs(g3_2_lin, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g4a = %f \n', eval(subs(g4a_lin, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g4b = %f \n', eval(subs(g4b_lin, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g5a = %f \n', eval(subs(g5a_lin, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g5b = %f \n', eval(subs(g5b_lin, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));


