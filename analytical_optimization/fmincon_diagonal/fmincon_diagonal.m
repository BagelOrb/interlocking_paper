clear all;
constants;
design_variables;

f =  @(x)  (wa(x) + wb(x)) *2* h / F(x);
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

problem.options = options;
problem.solver = 'fmincon';
problem.objective = f;
problem.nonlcon = @constraints;
problem.x0 = [1,1,1,1];

[x, fval] = fmincon(problem);
[c,ceq] = constraints(x);
 
fprintf('\n The minimum objective of %f with a max nominal stress of %f is reached for: \n', fval, 1 / fval)
fprintf('w_a* = %f \n', wa(x));
fprintf('w_b* = %f \n', wb(x));
fprintf('L* = %f \n', L(x));
fprintf('F* = %f \n', F(x));

fprintf('And occurs for the following constraint values: \n');
fprintf('g1a = %f \n', c(1));
fprintf('g1b = %f \n', c(2));
fprintf('g2 = %f \n', c(3));
fprintf('g3_1 = %f \n', c(4));
fprintf('g3_2 = %f \n', c(5));
fprintf('g4a = %f \n', c(6));
fprintf('g4b = %f \n', c(7));
fprintf('g5a = %f \n', c(8));
fprintf('g5b = %f \n', c(9));
fprintf('g6 = %f \n', c(10));
