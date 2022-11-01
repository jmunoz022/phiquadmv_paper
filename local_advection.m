function g = local_advection(xcoord,ycoord,beta,nquad)
%
% LOCAL_ADVECTION - Function that returns local contribution of the
%                   2D advection matrix
%
% INPUT: 
%   xcoord - physical coordinates in the x direction
%   ycoord - physical coordinates in the y direction
%   beta   - advection vector
%   nquad  - number of quadrature points
%
% OUTPUT: 
%   g - local contribution of the advection matrix
%

%Size of the element in each direction
a = xcoord(2)-xcoord(1);
b = ycoord(4)-ycoord(1);

%Quadrature points and weights
[z, w] = lobpts(nquad);

%Shape functions
N1 = @(x,y) (1-x)*(1-y)/4;
N2 = @(x,y) (1+x)*(1-y)/4;
N3 = @(x,y) (1+x)*(1+y)/4;
N4 = @(x,y) (1-x)*(1+y)/4;

%Gradients of shape function
N1grad = @(x,y) [(-(1-y)/4)*(2/a) (-(1-x)/4)*(2/b)];
N2grad = @(x,y) [((1-y)/4)*(2/a) (-(1+x)/4)*(2/b)];
N3grad = @(x,y) [((1+y)/4)*(2/a) ((1+x)/4)*(2/b)];
N4grad = @(x,y) [(-(1+y)/4)*(2/a) ((1-x)/4)*(2/b)];

%Jacobian
J = a*b/4;

%Initialize local matrix
g = zeros(4,4);

%Loop through quadrature points 
for i = 1:nquad
    for j = 1:nquad
        
        %Evaluate gradientes of shape functions times the advection vector
        bn1 = beta*N1grad(z(i),z(j))';
        bn2 = beta*N2grad(z(i),z(j))';
        bn3 = beta*N3grad(z(i),z(j))';
        bn4 = beta*N4grad(z(i),z(j))';
        
        %Evaluate shape functions
        n1 = N1(z(i),z(j));
        n2 = N2(z(i),z(j));
        n3 = N3(z(i),z(j));
        n4 = N4(z(i),z(j));
        
        %Accumulate local advection matrix
        g = g+[bn1*n1' bn1*n2' bn1*n3' bn1*n4';...
             bn2*n1' bn2*n2' bn2*n3' bn2*n4';...
             bn3*n1' bn3*n2' bn3*n3' bn3*n4';...
             bn4*n1' bn4*n2' bn4*n3' bn4*n4']*w(i)*w(j)*abs(J);
    end
end

end