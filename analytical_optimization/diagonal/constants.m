% Material A = PLA
% Material B = PP
% Units: 
% - Force: [N]
% - Dimensions: [mm]
% - Stress: [N/mm^2] 

E_a = 2797;
E_b = 302;

sigma_a = 47;
sigma_b = 10.5;

tau_a = sigma_a / sqrt(3);
tau_b = sigma_b / sqrt(3);

h = 0.5; % 0.5 mm finger height = average of previous opt.

w_a_min = 0.3;
w_b_min = 0.3;
w_max = 12 * 0.3;
L_min = 6 * 0.3;
L_max = 12 * 0.3;
