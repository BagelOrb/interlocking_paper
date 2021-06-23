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
NF = 20;

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
                x = [wa_, wb_, L_, F_];
                constraints = ...
                    [ g1a(x)
                      g1b(x)
                      g2(x)
                      g3_1(x)
                      g3_2(x)
                      g4a(x)
                      g4b(x)
                      g5a(x)
                      g5b(x)];
                constraints_are_valid = true;
                for g = constraints
                    constraints_are_valid = constraints_are_valid & (g <= 0.0001);
                end
                if constraints_are_valid
                    obj = f(x);
                    if obj < obj_min
                        obj_min = obj;
                        best_x = x;
                    end 
                   
                end
            end
        end
    end
end

x = best_x;
fprintf('The minimum objective of %f with a max nominal stress of %f is reached for: \n', obj_min, 1 / obj_min)
fprintf('w_a* = %f \n', wa(x));
fprintf('w_b* = %f \n', wb(x));
fprintf('L* = %f \n', L(x));
fprintf('F* = %f \n', F(x));

