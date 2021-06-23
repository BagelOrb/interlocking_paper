function c5b = g5b(F, w_a, w_b, L)
    constants;
    sigma_t = F/(w_b+h);
    sigma_m = (3*F*(w_a + w_b)*((w_a + w_b)^2 + 4*L^2))/( 8*h*w_b^2* L^2);
    tau = (3*F*(w_a + w_b))/(4*w_b*h*L);
    c5b = sqrt(((sigma_t + sigma_m)^2 )/2 + 3*tau^2)/sigma_b - 1;
end

