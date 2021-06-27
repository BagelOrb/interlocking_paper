function [A, W, dfdx] = getAandWmatrix(f, h, x, lambda)
    % Define matrix A
    
    h1 = h(1);
    h2 = h(2);
    h3 = h(3);
    h4 = h(4);
    
    dh1dx1 = diff(h1, x(1));
    dh1dx2 = diff(h1, x(2));
    dh1dx3 = diff(h1, x(3));
    dh1dx4 = diff(h1, x(4));

    dh2dx1 = diff(h2, x(1));
    dh2dx2 = diff(h2, x(2));
    dh2dx3 = diff(h2, x(3));
    dh2dx4 = diff(h2, x(4));

    dh3dx1 = diff(h3, x(1));
    dh3dx2 = diff(h3, x(2));
    dh3dx3 = diff(h3, x(3));
    dh3dx4 = diff(h3, x(4));

    dh4dx1 = diff(h4, x(1));
    dh4dx2 = diff(h4, x(2));
    dh4dx3 = diff(h4, x(3));
    dh4dx4 = diff(h4, x(4));

    A = [dh1dx1 dh1dx2 dh1dx3 dh1dx4;
         dh2dx1 dh2dx2 dh2dx3 dh2dx4;
         dh3dx1 dh3dx2 dh3dx3 dh3dx4;
         dh4dx1 dh4dx2 dh4dx3 dh4dx4];
     
%    
% 
%     A_lin = sym(size(A));
%     
%     
%     for i = 1:size(A, 1)
%         for j = 1:size(A,2)
%             A_lin(i,j) = linearize(A(i,j), x_eval(j), x(j));
%         end
%     end
%     
%             
    dfdx1 = diff(f, x(1));
    dfdx2 = diff(f, x(2));
    dfdx3 = diff(f, x(3));
    dfdx4 = diff(f, x(4));

    dfdx = [dfdx1;
            dfdx2;
            dfdx3;
            dfdx4];
 

    dLdx1 = dfdx1 + lambda.' * [dh1dx1; dh2dx1; dh3dx1; dh4dx1];
    dLdx2 = dfdx2 + lambda.' * [dh1dx2; dh2dx2; dh3dx2; dh4dx2] ;
    dLdx3 = dfdx3 + lambda.' * [dh1dx3; dh2dx3; dh3dx3; dh4dx3] ;
    dLdx4 = dfdx4 + lambda.' * [dh1dx4; dh2dx4; dh3dx4; dh4dx4];

    ddLdx1x1 = diff(dLdx1, x(1));   
    ddLdx1x2 = diff(dLdx1, x(2));
    ddLdx1x3 = diff(dLdx1, x(3));   
    ddLdx1x4 = diff(dLdx1, x(4));

    ddLdx2x1 = diff(dLdx2, x(1));   
    ddLdx2x2 = diff(dLdx2, x(2));
    ddLdx2x3 = diff(dLdx2, x(3));   
    ddLdx2x4 = diff(dLdx2, x(4));

    ddLdx3x1 = diff(dLdx3, x(1));   
    ddLdx3x2 = diff(dLdx3, x(2));
    ddLdx3x3 = diff(dLdx3, x(3));   
    ddLdx3x4 = diff(dLdx3, x(4));

    ddLdx4x1 = diff(dLdx4, x(1));   
    ddLdx4x2 = diff(dLdx4, x(2));
    ddLdx4x3 = diff(dLdx4, x(3));   
    ddLdx4x4 = diff(dLdx4, x(4));

    W = [ddLdx1x1 ddLdx1x2 ddLdx1x3 ddLdx1x4;
         ddLdx2x1 ddLdx2x2 ddLdx2x3 ddLdx2x4;
         ddLdx3x1 ddLdx3x2 ddLdx3x3 ddLdx3x4;
         ddLdx4x1 ddLdx4x2 ddLdx4x3 ddLdx4x4];
    
     
end
    