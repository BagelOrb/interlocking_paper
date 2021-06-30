clear all;

fvals = [];
lmax_vals = 1.2:.15:4.8;
for lmax = lmax_vals
    straight_constants;
    design_variables;
    objective;

    options = optimoptions('fmincon','Algorithm','sqp');

    problem.options = options;
    problem.solver = 'fmincon';
    problem.objective = f;
    problem.nonlcon = @(x) constraints(x, lmax);
    problem.x0 = [1,1,1,1,1];

    [x, fval] = fmincon(problem);
    fvals = [fvals fval];
    
end

max_stresses = 1 ./ fvals;
[max_stresses; lmax_vals]
figure(2);
plot(lmax_vals, max_stresses, '-o');
xlabel('$$\overline{L}$$ (mm)','interpreter','latex');
ylabel('Optimal ultimate stress (MPa)','interpreter','latex');
ax.FontSize = 15;
ylim([0 8]);

