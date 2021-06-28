function fun_quad = taylor2(fun, wa_k, wb_k, L_k, F_k)

    syms wa wb L F Z
    x = [wa; wb; L; F];
    x_k = [wa_k; wb_k; L_k; F_k];
    
    d_wa = diff(fun, wa);
    d_wb = diff(fun, wb);
    d_L = diff(fun, L);
    d_F = diff(fun, F);
    
    fun_lin = fun(wa_k, wb_k, L_k, F_k) + ...
        (wa-wa_k) * eval(subs(d_wa, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]) )+ ...
        (wb-wb_k) * eval(subs(d_wb, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]) )+ ...
        (L - L_k) * eval(subs(d_L, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]) ) + ...
        (F - F_k) * eval(subs(d_F, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]));
    
    [A, H, dfdx] =  getAandWmatrix(fun, [Z; Z; Z; Z], x, [Z; Z; Z; Z]);  
    H_k = subs(H, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]);
    
    fun_quad = fun_lin + (x - x_k).' * H_k * (x - x_k);
end

