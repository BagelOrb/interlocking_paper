w_a_list = linspace(0.3, 0.9, 7);
w_b_list = [0.3 0.4 0.5 0.6 0.7 0.8 0.9];
L_list = [6 7 8 9 10 11 12]*0.3;
F_list = 1:0.25:10;
obj_min = 9999999;

for w_a = w_a_list
    for w_b = w_b_list
        for L = L_list
            for F = F_list
                c1a = g1a(w_a);
                c1b = g1b(w_b);
                c2 = g2(w_a, w_b);
                c3_1 = g3_1(L);
                c3_2 = g3_2(L);
                c4a = g4a(F, w_a, w_b, L);
                c4b = g4b(F, w_a, w_b, L);
                c5a = g5a(F, w_a, w_b, L);
                c5b = g5b(F, w_a, w_b, L);
                if c1a <= 0 && c1b <= 0 && c2 <=0 && c3_1 <=0 && c3_2 <= 0 && c4a <= 0 && c4b <= 0 && c5a <= 0 && c5b<= 0
                    obj = f(w_a, w_b, F);
                    if obj < obj_min
                        obj_min = obj;
                        w_a_opt = w_a;
                        w_b_opt = w_b;
                        L_opt = L;
                        F_opt = F;
                    end 
                   
                end
            end
        end
    end
end

fprintf('The minimum objective of %f is reached for: \n', obj_min)
fprintf('w_a* = %f \n', w_a_opt);
fprintf('w_b* = %f \n', w_b_opt);
fprintf('L* = %f \n', L_opt);
fprintf('F* = %f \n', F_opt);
            

