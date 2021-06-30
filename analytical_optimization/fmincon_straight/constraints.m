function [c,ceq] = constraints(x)
    straight_constants;
    design_variables;

    c(1) = 1 - wb(x) / 2 / wb_min;
    c(2) = 1 - va(x) / wa_min;
    c(3) = 1 - vb(x) / wb_min;
    c(4) = 1 - hf(x) / h_min;
    c(5) =  (va(x) + vb(x)) / L_max - 1;
    c(6) = 1 - wa * hf(x) * sa / F(x);
    c(7) = 1 - wb(x) * hf(x) * sb / F(x);
    c(8) = 1 - 2 * va(x) * sa / (F(x) * sqrt(3 * (wb(x)/(wa+wb(x))/hc)^2 + 3*(sa/saz/wa)^2));
    c(9) = 1 - 2*vb(x)*wb(x)*tbz / F(x);

    ceq = [];
end
