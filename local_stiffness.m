function k = local_stiffness(xcoord,ycoord,diff,nquad)
%
% LOCAL_STIFFNESS - Function that returns local contribution of the
%                   2D stiffness matrix
%
% INPUT: 
%   xcoord - physical coordinates in the x direction
%   ycoord - physical coordinates in the y direction
%   diff   - diffusion coeficient 
%   nquad  - number of quadrature points
%
% OUTPUT: 
%   k - local contribution of the stiffness matrix
%

%Size of the element in each direction
a = xcoord(2)-xcoord(1);
b = ycoord(4)-ycoord(1);

%Quadrature points and weights
[z, w] = lobpts(nquad);

%Gradient of shape functions
N1grad = @(x,y) [(-(1-y)/4)*(2/a) (-(1-x)/4)*(2/b)];
N2grad = @(x,y) [((1-y)/4)*(2/a) (-(1+x)/4)*(2/b)];
N3grad = @(x,y) [((1+y)/4)*(2/a) ((1+x)/4)*(2/b)];
N4grad = @(x,y) [(-(1+y)/4)*(2/a) ((1-x)/4)*(2/b)];

%%%Jacobian%%%
J = a*b/4;

%Initialize local matrix
k = zeros(4,4);

%Loop through quadrature points 
for i = 1:nquad
    for j = 1:nquad
        
        %Evaluate gradientes of shape functions
        n1 = N1grad(z(i),z(j));
        n2 = N2grad(z(i),z(j));
        n3 = N3grad(z(i),z(j));
        n4 = N4grad(z(i),z(j));
        
        %Accumulate local stiffness matrix
        k = k+[n1*n1' n1*n2' n1*n3' n1*n4';...
             n2*n1' n2*n2' n2*n3' n2*n4';...
             n3*n1' n3*n2' n3*n3' n3*n4';...
             n4*n1' n4*n2' n4*n3' n4*n4']*diff*w(i)*w(j)*abs(J);
    end
end

end