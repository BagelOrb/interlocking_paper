clear all;
straight_case;

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

dgdx_latex = string([]);
for j = 1:size(dgdx,2)
    for i = 1:size(dgdx,1)
        p = dgdx(i,j);
        dgdx_latex(i,j) = string(latex(p));
    end
end

dfdx_latex = string([]);
for i = 1:length(dfdx)
    p = dfdx(i);
    dfdx_latex(i) = string(latex(p));
end

dLgdLx_latex = string([]);
for j = 1:size(dLgdLx,2)
    for i = 1:size(dLgdLx,1)
        p = dLgdLx(i,j);
        dLgdLx_latex(i,j) = string(latex(p));
    end
end

dLfdLx_latex = string([]);
for i = 1:length(dLfdLx)
    p = dLfdLx(i);
    dLfdLx_latex(i) = string(latex(p));
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