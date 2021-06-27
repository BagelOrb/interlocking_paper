clear constants; % otherwise changes to the script aren't loaded until you restart MATLAB
constants;
clear design_variables;
design_variables;
clear objective;
objective;
clear constraints;
constraints;

Nwa = 10;
Nwb = 10;
NF = 10;

was = linspace(2 * wa_min, w_max - 2 * wa_min, Nwa);
wbs = linspace(2 * wb_min, w_max - 2 * wb_min, Nwb);
Fs = linspace(1, 10, NF);
Ls = L_min:0.3:L_max; 

obj_min = 9999999;
best_x = [-1,-1,-1,-1]; % placeholder


for wa_ = was
    for wb_ = wbs
        for L_ = Ls
            for F_ = Fs
                x_k = [wa_, wb_, L_, F_];
              
                constraints = ...
                    [ g1a(x_k)
                      g1b(x_k)
                      g2(x_k)
                      g3_1(x_k)
                      g3_2(x_k)
                      g4a(x_k)
                      g4b(x_k)
                      g5a(x_k)
                      g5b(x)];
                  
         
                constraints_are_valid = true;
                for g = constraints
                    constraints_are_valid = constraints_are_valid & (g <= 0.0001);
                end
                if constraints_are_valid
                    obj = f(x_k);
                    if obj < obj_min
                        obj_min = obj;
                        best_x = x_k;
                    end 
                   
                end
            end
        end
    end
end

x_k = best_x;
fprintf('The minimum objective of %f with a max nominal stress of %f is reached for: \n', obj_min, 1 / obj_min)
fprintf('w_a* = %f \n', wa(x_k));
fprintf('w_b* = %f \n', wb(x_k));
fprintf('L* = %f \n', L(x_k));
fprintf('F* = %f \n', F(x_k));


fprintf('And occurs for the following constraint values: \n');
fprintf('g1a = %f \n', g1a(x_k));
fprintf('g1b = %f \n', g1b(x_k));
fprintf('g2 = %f \n', g2(x_k));
fprintf('g3_1 = %f \n', g3_1(x_k));
fprintf('g3_2 = %f \n', g3_2(x_k));
fprintf('g4a = %f \n', g4a(x_k));
fprintf('g4b = %f \n', g4b(x_k));
fprintf('g5a = %f \n', g5a(x_k));
fprintf('g5b = %f \n', g5b(x_k));



