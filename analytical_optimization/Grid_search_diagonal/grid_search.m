clear constants; % otherwise changes to the script aren't loaded until you restart MATLAB
constants;
clear design_variables;
design_variables;
clear objective;
objective;
clear constraints;
constraints;

syms wa wb L F

g1a = sym(g1a);
g1b = sym(g1b);
g2 = sym(g2);
g3_1 = sym(g3_1);
g3_2 = sym(g3_2);
g4a = sym(g4a);
g4b = sym(g4b);
g5a = sym(g5a);
g5b = sym(g5b);
g6 = sym(g6);


if exist('iter','var') == 0
    Niter = 10;
    was = linspace(2 * wa_min, w_max - 2 * wa_min, Niter);
    wbs = linspace(2 * wb_min, w_max - 2 * wb_min, Niter);
    Fs = linspace(1, 10, Niter);
    Ls = L_min:0.3:L_max; 

else
    Niter = 5;
    was = linspace(wa_best - 0.3/(iter^2), wa_best + 0.3/(iter^2), Niter);
    wbs = linspace(wb_best - 0.3/(iter^2), wb_best + 0.3/(iter^2), Niter);
    Fs = linspace(F_best - 0.5/(iter^2), F_best + 0.5/(iter^2), Niter);
    Ls = linspace(L_best - 0.5/(iter^2), L_best + 0.5/(iter^2), Niter);


end

obj_min = 9999999;
best_x = [-1,-1,-1,-1]; % placeholder
i = 0;

for wa = was
    i = i+1;
    disp(['Start iteration ', num2str(i), ' out of ', num2str(Niter)]);
    for wb = wbs
        for L = Ls
            for F = Fs
                x = [wa, wb, L, F];
                
                
                constraints = [ ...
                    double(eval(g1a)) ...
                    double(eval(g1b)) ...
                    double(eval(g2)) ...
                    double(eval(g3_1)) ...
                    double(eval(g3_2)) ...
                    double(eval(g4a)) ...
                    double(eval(g4b)) ...
                    double(eval(g5a)) ...
                    double(eval(g5b)) ... 
                    double(eval(g6)) ...
                     ];
                  
              
                    
                constraints_are_valid = true;
                
                for g = constraints
                    constraints_are_valid = constraints_are_valid & (g <= 0.00000001);
                end
                
                if constraints_are_valid
                    obj_eval = double(eval(subs(f)));
                    if obj_eval < obj_min
                        obj_min = obj_eval;
                        wa_best = wa;
                        wb_best = wb;
                        L_best = L;
                        F_best = F;
                      
                    end 
                   
                end
            end
        end
    end
end

wa = wa_best;
wb = wb_best;
L = L_best;
F = F_best;

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
fprintf('g1a = %f \n', eval(subs(g1a)));
fprintf('g1b = %f \n', eval(subs(g1b)));
fprintf('g2 = %f \n', eval(subs(g2)));
fprintf('g3_1 = %f \n', eval(subs(g3_1)));
fprintf('g3_2 = %f \n', eval(subs(g3_2)));
fprintf('g4a = %f \n', eval(subs(g4a)));
fprintf('g4b = %f \n', eval(subs(g4b)));
fprintf('g5a = %f \n', eval(subs(g5a)));
fprintf('g5b = %f \n', eval(subs(g5b)));
fprintf('g6 = %f \n', eval(subs(g6)));



