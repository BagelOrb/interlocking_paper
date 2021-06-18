constants;
design_variables;

gwb = @(x) 1 - wb(x) / 2 / wb_min;
gva = @(x) 1 - va(x) / wa_min;
gvb = @(x) 1 - vb(x) / wb_min;
ghf = @(x) 1 - hf(x) / h_min;
gd =  @(x) (va(x) + vb(x)) / L_max - 1;
gta = @(x) 1 - wa * hf(x) * sa / F(x);
gtb = @(x) 1 - wb(x) * hf(x) * sb / F(x);
gca = @(x) 1 - 2 * va(x) * sa / (F(x) * sqrt(3 * (wb(x)/(wa+wb(x))/hc)^2 + 3*(sa/saz/wa)^2));
gzb = @(x) 1 - 2*vb(x)*wb(x)*tbz / F(x);
