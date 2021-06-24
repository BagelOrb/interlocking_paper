function fun_lin = linearize(fun, x_eval, x)
    fun_lin = subs(fun, x_eval) +   (x-x_eval) * subs(diff(fun),x_eval);
end

