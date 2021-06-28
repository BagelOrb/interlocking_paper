function [A, W, dfdx] = getAandWmatrix(f, h, x, lambda)

    nx = size(x,1);
    nl = size(lambda,1);
    lagrange = sym(f);
    for i = 1:nl
        lagrange = lagrange + lambda(i) * h(i);
    end
    symbols = [x; lambda];
    n = size(symbols,1);
    mat = sym('a',[nx, n]);
    dfds = [];
    for x = 1:n
        dfds = [dfds ; diff(lagrange, symbols(x))];
    end
    for x = 1:nx
        for y = 1:n
            mat(y,x) = diff(dfds(x), symbols(y));
        end
    end
    
    W = mat(1:nx,1:nx);
    A = mat(nx+1:nx+nl,1:nx);
    dfdx = dfds(1:nx);
end
    