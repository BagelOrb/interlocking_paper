clear; % otherwise changes to the script aren't loaded until you restart MATLAB
constants;
design_variables;
objective_sqp;
constraints_sqp;

syms wa wb L F
x_eval = [1; 1; 1; 1];
x = [wa; wb; L; F];
obj_min = 999999;

opts = optimoptions(@quadprog,'Algorithm','trust-region-reflective');
wa = x_eval(1);
wb = x_eval(2);
L = x_eval(3);
F = x_eval(4);

n = length(x);
n_iter = 100;
delta = 10000;
lambda = 0.01*ones(n,1);

% Set active constraints
h1 = g3_2;
h2 = g5a;
h3 = g5b;
h4 = g6;
h = [sym(h1); sym(h2); sym(h3); sym(h4)];

for n = 1:n_iter
    
    % lb = [wa_min - x_eval(1) , wb_min -  x_eval(2), L_min - x_eval(3), 99999]
    % ub = [w_max - x_eval(1), w_max - x_eval(2), L_max, -99999]

    lambda = double(lambda + 1/delta*subs(h));
   
    [A, W, dfdx] =  getAandWmatrix(f, h, x, lambda);
    
    d = x - x_eval;
    W_eval = double(subs(W));
    dfdx_eval = double(subs(dfdx));
    A_eval = double(subs(A));
    h_eval = double(subs(h));
    
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
    x_eval = x_eval + dx;

    ev = eig(W_eval)
    isposdef = all(ev>-0);
    
   
    wa = x_eval(1);
    wb = x_eval(2);
    L = x_eval(3);
    F = x_eval(4);
    
     if isposdef
        obj_min  = f(wa, wb, L, F)
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



