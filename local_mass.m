function m = local_mass(xcoord,ycoord,nquad)
%
% LOCAL_MASS - Function that returns local contribution of the
%              2D mass matrix
%
% INPUT: 
%   xcoord - physical coordinates in the x direction
%   ycoord - physical coordinates in the y direction
%   nquad  - number of quadrature points
%
% OUTPUT: 
%   m - local contribution of the mass matrix
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

%Jacobian
J = a*b/4;

%Initialize local matrix
m = zeros(4,4);

%Loop through quadrature points 
for i = 1:nquad
    for j = 1:nquad
        
        %Evaluate shape functions
        n = [N1(z(i),z(j)) N2(z(i),z(j)) N3(z(i),z(j)) N4(z(i),z(j))];
        
        %Accumulate local mass matrix
        m = m+n'*n*w(i)*w(j)*abs(J);
    end
end

end
