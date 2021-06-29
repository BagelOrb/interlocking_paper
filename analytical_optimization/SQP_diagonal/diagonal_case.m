constants;

syms wa wb L F;
x = [wa; wb; L; F];

f = (wa + wb) *2* h  / F;

g1a = 1- wa/wa_min;
g1b = 1- wb/wb_min;
g2  =(wa + wb)/w_max -1;
g3_1 = 1- L/L_min;
g3_2 = L/L_max - 1;
g4a = (3*F*(wa+wb))/(taz*wa^2*L) -1 ;
g4b = (3*F*(wa+wb))/(tbz*wb^2*L) -1 ;
g5a = constraint5(F, wa, wb, L, sa);
g5b = constraint5(F, wb, wa, L, sb);
g6 = 1- F/F_min;

g = [ g1a, g1b, g2, g3_1, g3_2, g4a, g4b, g5a, g5b, g6 ];
g_names = [ "g1a", "g1b", "g2", "g3_1", "g3_2", "g4a", "g4b", "g5a", "g5b", "g6" ];

function g5 = constraint5(F, w1, w2, L, s1)
    constants;
    sigma_t = F/(w1+h);
    sigma_m = (3*F*(w1 + w2)*((w1 + w2)^2 + 4*L^2))/( 8*h*w1^2* L^2);
    tau = (3*F*(w1 + w2))/(4*w1*h*L);
    g5 = sqrt(((sigma_t + sigma_m)^2 )/2 + 3*tau^2)/s1 - 1;
end