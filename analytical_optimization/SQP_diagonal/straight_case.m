straight_constants;

syms wb va vb hf F;
x = [wb va vb hf F].';

f = (wa + wb) * (hf + hc) / F;

gwb = 1 - wb / 2 / wb_min;
gva = 1 - va / wa_min;
gvb = 1 - vb / wb_min;
ghf = 1 - hf / h_min;
gd =  (va + vb) / L_max - 1;
gta = 1 - wa * hf * sa / F;
gtb = 1 - wb * hf * sb / F;
gca = 1 - 2 * va * sa / (F * sqrt(3 * (wb/(wa+wb)/hc)^2 + 3*(sa/saz/wa)^2));
gzb = 1 - 2*vb*wb*tbz / F;

g = [ gwb gva gvb ghf gd gta gtb gca gzb ];
g_names = [ "gwb" "gva" "gvb" "ghf" "gd" "gta" "gtb" "gca" "gzb" ];