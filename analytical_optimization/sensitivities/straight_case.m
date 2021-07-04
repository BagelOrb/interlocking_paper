
move_limit = 1; 

syms wb va vb hf F;
syms hc wa wa_min wb_min sa sb sbz saz tbz L_max h_min 

x = [wb; va; vb; hf; F];
nx = length(x);
%x_k = [2.7, 2.16, 1.44, 1.0, 28.2]; % starting point
x_k = ones(1,nx); % starting point

f = (wa + wb) * (hf + hc) / F;

gwb = 1 - wb / 2 / wb_min;
gva = 1 - va / wa_min;
gvb = 1 - vb / wb_min;
ghf = 1 - hf / h_min;
gd =  (va + vb) / L_max - 1;
gta = 1 - wa * hf * sa / F;
gtb = 1 - wb * hf * sb / F;
gca = 1 - 2 * va * sa * (wa+wb) / (F * sqrt(3 * (wb/hc)^2 + 3*(sa/saz)^2));
gzb = 1 - 2*vb*wb*tbz / F;

g = [ gwb, gva, gvb, ghf, gd, gta, gtb, gca, gzb ];
g_names = [ "g1", "g2a", "g2b", "g3", "4", "g5a", "g5b", "g6", "g7" ];







