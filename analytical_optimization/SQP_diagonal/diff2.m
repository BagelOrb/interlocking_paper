function [ddfdss, dfds] = diff2(fun, symbols)
    n = size(symbols, 1);
    ddfdss = sym('a',n);
    dfds = [];
    for x = 1:n
        dfds = [dfds ; diff(fun, symbols(x))];
    end
    for x = 1:n
        for y = 1:n
            ddfdss(y, x) = diff(dfds(x), symbols(y));
        end
    end
end