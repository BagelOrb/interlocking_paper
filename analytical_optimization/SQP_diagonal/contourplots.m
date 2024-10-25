function contourplots(x, x_history, x1_array, x2_array, f_f, g_f, g_names)
% Initialization
clf, hold off
warning('off')
% Combinations of design variables D and d 
nx = length(x);
n_constr = length(g_f(1, 1));
g_eval = zeros(length(x2_array), length(x1_array), n_constr);
% Matrix of output values for combinations of design variables D and d: 


if nx > 2
    for r = 3:nx
        x_k(r) = x_history(size(x_history, 1),r);
    end
end

for j=1:1:length(x2_array)
  for i=1:1:length(x1_array)
 
    % Assignment of design variables:
    x_k(1) = x1_array(i);
    x_k(2) = x2_array(j);
    
    x_kc = num2cell(x_k, 1);
    
    % x_k = [x_k(1) x_k(2)];
    
 	 % Objective function
    f_eval_k = f_f(x_kc{:});
    g_eval_k = g_f(x_kc{:});
    
    
    if f_eval_k > 10
        f_eval_k = (f_eval_k)/10 + 9;    % Set maximum to better visualize levels in contour plot
    end
    
    % Grid value of objective function:
    fobj(j,i) = f_eval_k; 
    
    set(gcf, 'WindowState', 'maximized');
    
    for p = 1:n_constr
    % Grid values of constraints:
      g_eval(j,i,p) = g_eval_k(p);    % Scaled length constraint
      
    
    end
  end
end

% Contour plot of scaled spring problem
contourf(x1_array, x2_array, fobj, 30,'LineColor','none');
colormap('parula')
% contour(x1_array, x2_array, fobj, 'ShowText','on');
xlabel(string(x(1))), ylabel(string(x(2))), ...
   %title('Contour plot')
hold on

for p = 1:n_constr
    [C, hContour] = contour(x1_array, x2_array, g_eval(:,:,p), [0.0 0.0], 'LineColor', [1 1 1]*230/255,  'showtext', 'on');
    clabel(C,hContour,'Color', [1 1 1]*230/255)
    drawnow;
    lab = char(g_names(p));
    labels=hContour.TextPrims;
    for idx = 1 : numel(labels)
        hContour.TextPrims(idx).String= lab;
    end
    contour(x1_array, x2_array, g_eval(:,:,p), [0.02, 0.02],'--', 'LineColor', [1 1 1]*230/255) ;% Infeasible region
    
    hold on
end

plot(x_history(:,1), x_history(:,2), 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 2);
plot(x_history(size(x_history, 1),1), x_history(size(x_history, 1),2), 'Marker', 'p', 'MarkerFaceColor', 'r','MarkerEdgeColor', 'b', 'MarkerSize', 20);
grid

end 