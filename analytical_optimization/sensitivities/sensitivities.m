clear all;
diagonal_case;

ng = length(g);
nx = length(x);

dgdx = sym(zeros(nx, ng));
dLgdLx = sym(zeros(nx, ng));

[dLfdLx, dfdx] =  sens(f, x);

for i = 1:ng
    [dLgdLx_i, dgdx_i] = sens(g(i), x);
    dgdx(:,i) = dgdx_i;
    dLgdLx(:,i) = dLgdLx_i;
end



function [dLfdLx, dfdx] = sens(f, x)
    n = size(x, 1);
    dfdx = [];
    dLfdLx = [];
    
    for i = 1:n
        dfdx = [dfdx ; simplify(diff(f, x(i)))];
    end
    for i = 1:n
        dLfdLx = [dLfdLx; simplify(x(i) / f * dfdx(i))];
    end
end