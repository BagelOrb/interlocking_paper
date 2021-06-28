function [c,ceq] = constraints(x)
    constants;
    design_variables;
    c(1) =  1- wa(x)/wa_min; 
    c(2) =  1- wb(x)/wb_min; 
    c(3)  = (wa(x) + wb(x))/w_max -1; 
    c(4) =  1- L(x)/L_min; 
    c(5) = L(x)/L_max - 1; 
    c(6) = (3*F(x)*(wa(x)+wb(x)))/(taz*wa(x)^2*L(x)) -1 ; 
    c(7) = (3*F(x)*(wa(x)+wb(x)))/(tbz*wb(x)^2*L(x)) -1 ; 
    c(8) = constraint5(F(x), wa(x), wb(x), L(x), sa); 
    c(9) = constraint5(F(x), wb(x), wa(x), L(x), sb); 
    c(10) = 1- F(x)/F_min; 
    
    ceq = [];
end


function g5 = constraint5(F, w1, w2, L, s1)
    constants;
    sigma_t = F/(w1+h);
    sigma_m = (3*F*(w1 + w2)*((w1 + w2)^2 + 4*L^2))/( 8*h*w1^2* L^2);
    tau = (3*F*(w1 + w2))/(4*w1*h*L);
    g5 = sqrt(((sigma_t + sigma_m)^2 )/2 + 3*tau^2)/s1 - 1;
end