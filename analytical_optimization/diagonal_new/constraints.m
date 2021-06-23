constants;
design_variables;

g1a = @(x) 1- wa(x)/wa_min;
g1b = @(x) 1- wb(x)/wb_min;
g2  = @(x)(wa(x) + wb(x))/w_max -1;
g3_1 = @(x) 1- L(x)/L_min;
g3_2 = @(x) L(x)/L_max - 1;
g4a = @(x) (3*F(x)*(wa(x)+wb(x)))/(taz*wa(x)^2*L(x)) -1 ;
g4b = @(x) (3*F(x)*(wa(x)+wb(x)))/(tbz*wb(x)^2*L(x)) -1 ;
g5a = @(x) constraint5(F(x), wa(x), wb(x), L(x), sa);
g5b = @(x) constraint5(F(x), wb(x), wa(x), L(x), sb);


function g5 = constraint5(F, w1, w2, L, s1)
    constants;
    sigma_t = F/(w1+h);
    sigma_m = (3*F*(w1 + w2)*((w1 + w2)^2 + 4*L^2))/( 8*h*w1^2* L^2);
    tau = (3*F*(w1 + w2))/(4*w1*h*L);
    g5 = sqrt(((sigma_t + sigma_m)^2 )/2 + 3*tau^2)/s1 - 1;
end