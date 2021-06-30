clear all;
straight_constants;
design_variables;
objective;

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

problem.options = options;
problem.solver = 'fmincon';
problem.objective = f;
problem.nonlcon = @constraints;
problem.x0 = [1,1,1,1,1];

[x, fval] = fmincon(problem);
[c,ceq] = constraints(x);
 
fprintf('\n The minimum objective of %f with a max nominal stress of %f is reached for: \n', fval, 1 / fval)
fprintf('wb = %f \n', wb(x));
fprintf('va = %f \n', va(x));
fprintf('vb = %f \n', vb(x));
fprintf('hf = %f \n', hf(x));
fprintf('F* = %f \n', F(x));

fprintf('And occurs for the following constraint values: \n');
fprintf('c1 = %f \n', c(1));
fprintf('c2 = %f \n', c(2));
fprintf('c3 = %f \n', c(3));
fprintf('c4 = %f \n', c(4));
fprintf('c5 = %f \n', c(5));
fprintf('c6 = %f \n', c(6));
fprintf('c7 = %f \n', c(7));
fprintf('c8 = %f \n', c(8));
fprintf('c9 = %f \n', c(9));
