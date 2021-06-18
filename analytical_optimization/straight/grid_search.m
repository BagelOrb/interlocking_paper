clear constants; % otherwise changes to the script aren't loaded until you restart MATLAB
constants;
clear design_variables;
design_variables;
clear objective;
objective;
clear constraints;
constraints;

Nwb = 45;
Nva = 45;
Nhf = 13;
NF = 50;

wbs = linspace(2 * wb_min, w_max - 2 * wa_min, Nwb);
vas = linspace(wa_min, w_max - wb_min, Nva);
hfs = linspace(h_min, h_max - h_min, Nhf);
Fs = linspace(10, 30, NF);

obj_min = 9999999;
best_x = [-1,-1,-1,-1]; % placeholder

for wb_ = wbs
    for va_ = vas
        for hf_ = hfs
            for F_ = Fs
                x = [wb_, va_, hf_, F_];
                constraints = ...
                    [ gwb(x)
                      gva(x)
                      gvb(x)
                      ghf(x)
                      gd(x)
                      gta(x)
                      gtb(x)
                      gca(x)
                      gzb(x)];
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
fprintf('w_b* = %f \n', wb(x));
fprintf('v_a* = %f \n', va(x));
fprintf('h_f* = %f \n', hf(x));
fprintf('F* = %f \n', F(x));
fprintf('v_b* = %f \n', vb(x));
fprintf('w_a* = %f \n', wa);
fprintf('h_c* = %f \n', hc);
            

