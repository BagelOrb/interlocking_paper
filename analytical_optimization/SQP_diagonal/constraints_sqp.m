constants;
design_variables;

g1a = @(wa, wb, L, F) 1- wa/wa_min;
g1b = @(wa, wb, L, F) 1- wb/wb_min;
g2  = @(wa, wb, L, F)(wa + wb)/w_max -1;
g3_1 = @(wa, wb, L, F) 1- L/L_min;
g3_2 = @(wa, wb, L, F) L/L_max - 1;
g4a = @(wa, wb, L, F) (3*F*(wa+wb))/(taz*wa^2*L) -1 ;
g4b = @(wa, wb, L, F) (3*F*(wa+wb))/(tbz*wb^2*L) -1 ;
g5a = @(wa, wb, L, F) constraint5(F, wa, wb, L, sa);
g5b = @(wa, wb, L, F) constraint5(F, wb, wa, L, sb);
g6 = @(wa, wb, L, F) 1- F/F_min;

function g5 = constraint5(F, w1, w2, L, s1)
    constants;
    sigma_t = F/(w1+h);
    sigma_m = (3*F*(w1 + w2)*((w1 + w2)^2 + 4*L^2))/( 8*h*w1^2* L^2);
    tau = (3*F*(w1 + w2))/(4*w1*h*L);
    g5 = sqrt(((sigma_t + sigma_m)^2 )/2 + 3*tau^2)/s1 - 1;
end