clear constants; % otherwise changes to the script aren't loaded until you restart MATLAB
constants;
clear design_variables;
design_variables;
clear objective;
objective;
clear constraints;
constraints;

syms wa wb L F
syms wa_k wb_k L_k F_k

if exist('iter','var') == 0
    Niter = 10;
    was = linspace(2 * wa_min, w_max - 2 * wa_min, Niter);
    wbs = linspace(2 * wb_min, w_max - 2 * wb_min, Niter);
    Fs = linspace(1, 5, Niter);
    Ls = L_min:0.3:L_max; 

else
    Niter = 10;
    was = linspace(wa_best - 0.3/(iter^2), wa_best + 0.3/(iter^2), Niter);
    wbs = linspace(wb_best - 0.3/(iter^2), wb_best + 0.3/(iter^2), Niter);
    Fs = linspace(F_best - 0.5/(iter^2), F_best + 0.5/(iter^2), Niter);
    Ls = linspace(L_best - 0.5/(iter^2), L_best + 0.5/(iter^2), Niter);


end

obj_min = 9999999;
best_x = [-1,-1,-1,-1]; % placeholder
i = 0;


disp('Approximating Taylor series...');

g1a_lin =  taylor2(g1a, wa_k, wb_k, L_k, F_k);
fprintf('=');pause(0.001);
g1b_lin = taylor2(g1b, wa_k, wb_k, L_k, F_k);
fprintf('=');pause(0.001);
g2_lin = taylor2(g2, wa_k, wb_k, L_k, F_k);
fprintf('=');pause(0.001);
g3_1_lin = taylor2(g3_1, wa_k, wb_k, L_k, F_k);
fprintf('=');pause(0.001);
g3_2_lin = taylor2(g3_2, wa_k, wb_k, L_k, F_k);
fprintf('=');pause(0.001);
g4a_lin = taylor2(g4a, wa_k, wb_k, L_k, F_k);
fprintf('=');pause(0.001);
g4b_lin = taylor2(g4b, wa_k, wb_k, L_k, F_k);
fprintf('=');pause(0.001);
g5a_lin = taylor2(g5a, wa_k, wb_k, L_k, F_k);
fprintf('=');pause(0.001);
g5b_lin = taylor2(g5b, wa_k, wb_k, L_k, F_k);
fprintf('=');pause(0.001);
f_lin = taylor2(f, wa_k, wb_k, L_k, F_k);
fprintf('= \n');

g1a_k = subs(g1a_lin, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]);
g1b_k = subs(g1b_lin, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]);
g2_k = subs(g2_lin, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]);
g3_1_k = subs(g3_1_lin, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]);
g3_2_k = subs(g3_2_lin, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]);
g4a_k = subs(g4a_lin, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]);
g4b_k = subs(g4b_lin, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]);
g5a_k = subs(g5a_lin, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]);
g5b_k = subs(g5b_lin, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]);
f_k = subs(f_lin, [wa, wb, L, F], [wa_k, wb_k, L_k, F_k]);



for wa_k = was
    i = i+1;
    disp(['Start iteration ', num2str(i), ' out of ', num2str(Niter)]);
    for wb_k = wbs
        for L_k = Ls
            for F_k = Fs
                x_k = [wa_k, wb_k, L_k, F_k];
                
                
                constraints = [ ...
                    double(eval(g1a_k)) ...
                    double(eval(g1b_k)) ...
                    double(eval(g2_k)) ...
                    double(eval(g3_1_k)) ...
                    double(eval(g3_2_k)) ...
                    double(eval(g4a_k)) ...
                    double(eval(g4b_k)) ...
                    double(eval(g5a_k)) ...
                    double(eval(g5b_k)) ...              
                     ];
                  
              
                    
                constraints_are_valid = true;
                
                for g = constraints
                    constraints_are_valid = constraints_are_valid & (g <= 0.00000001);
                end
                
                
                
                
                if constraints_are_valid
                    obj_eval = double(eval(f_k));
                    if obj_eval < obj_min
                        obj_min = obj_eval;
                        wa_best = wa_k;
                        wb_best = wb_k;
                        L_best = L_k;
                        F_best = F_k;
                      
                    end 
                   
                end
            end
        end
    end
end


if exist('iter','var') == 1
    iter = iter + 1
else 
    iter = 1
    
end 

x_k = best_x;
fprintf('The minimum objective of %f with a max nominal stress of %f is reached for: \n', obj_min, 1 / obj_min)
fprintf('w_a* = %f \n', wa_best);
fprintf('w_b* = %f \n', wb_best);
fprintf('L* = %f \n', L_best);
fprintf('F* = %f \n', F_best);


fprintf('And occurs for the following constraint values: \n');
fprintf('g1a = %f \n', eval(subs(g1a_lin, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g1b = %f \n', eval(subs(g1b_lin, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g2 = %f \n', eval(subs(g2_lin, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g3_1 = %f \n', eval(subs(g3_1_lin, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g3_2 = %f \n', eval(subs(g3_2_lin, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g4a = %f \n', eval(subs(g4a_lin, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g4b = %f \n', eval(subs(g4b_lin, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g5a = %f \n', eval(subs(g5a_lin, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));
fprintf('g5b = %f \n', eval(subs(g5b_lin, [wa, wb, L, F], [wa_best, wb_best, L_best, F_best])));



