clear; % otherwise changes to the script aren't loaded until you restart MATLAB
constants;
design_variables;
objective_sqp;
constraints_sqp;

syms wa wb L F
x_eval = [0.6; 1.4; 2.7; 3.0];
x = [wa; wb; L; F];

opts = optimoptions(@quadprog,'Algorithm','active-set');
wa = x_eval(1);
wb = x_eval(2);
L = x_eval(3);
F = x_eval(4);

n = length(x);
n_iter = 10;
lambda = ones(n,1);

% Set active constraints
h1 = g3_2;
h2 = g4a;
h3 = g5a;
h4 = g5b;
h = [sym(h1); sym(h2); sym(h3); sym(h4)];

for n = 1:n_iter
    
    % lb = [wa_min - x_eval(1) , wb_min -  x_eval(2), L_min - x_eval(3), 99999]
    % ub = [w_max - x_eval(1), w_max - x_eval(2), L_max, -99999]

    
    [A, W, dfdx] =  getAandWmatrix(f, h1, h2, h3, h4, x, x_eval, lambda);

    W_eval = double(subs(W));
    dfdx_eval = double(subs(dfdx));
    A_eval = double(subs(A));
    h_eval = double(subs(h));

    dx = quadprog(W_eval, dfdx_eval.',[], [],  A_eval, -h_eval);
    lambda = rref(dx.'*A_eval.');

    x_eval = x_eval + dx;

    d = eig(W_eval);
    isposdef = all(d>0);
    
    
    wa = x_eval(1);
    wb = x_eval(2);
    L = x_eval(3);
    F = x_eval(4);
    
    obj_min  = f(wa, wb, L, F)

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



