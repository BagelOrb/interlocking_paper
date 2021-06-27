% Material A = PLA
% Material B = PP
% Units: 
% - Force: [N]
% - Dimensions: [mm]
% - Stress: [N/mm^2] 

E_a = 2797;
E_b = 302;

sa = 47;
sb = 10.5;
saz = 33;
sbz = 10.6;

ta = sa / sqrt(3);
tb = sb / sqrt(3);
taz = saz / sqrt(3);
tbz = sbz / sqrt(3);

h = 0.5; % 0.5 mm finger height = average of previous opt.

h_min = 0.2;
h_max = 12 * h_min;

wa_min = 0.3;
wb_min = 0.3;
w_max = 12 * 0.3;
L_max = 12 * 0.3;
L_min = 6 * 0.3;
F_min = 1;
