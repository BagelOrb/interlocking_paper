function fun_lin = linearize(fun, x_k, x)
    fun_lin = subs(fun, x_k) +   (x-x_k) * subs(diff(fun),x_k);
end

