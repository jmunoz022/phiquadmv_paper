function dir = BC(nx,ny,P)
%
% BC - Function that returns Dirichlet nodes for each problem
%                   
% INPUT: 
%   nx - number of nodes in the x direction
%   ny - number of nodes in the x direction
%   P  - problem number:
%        - P=1: Heat equation
%        - P=2: Eriksson-Johnson problem 
%        - P=3: Hochbruck-Osterman equation
%        - P=4: Allen-Cahn equation 
%
% OUTPUT: 
%   dir - vector containing the Dirichlet nodes
%

if P == 1 %Heat equation
    dir = zeros(1,2*ny+2*nx-4);
    dir(1:nx) = 1:nx;
    for i = 1:ny-2
       dir(nx+2*i-1) = i*nx+1;
       dir(nx+2*i) = (i+1)*nx;
    end
    dir(nx+2*(ny-2)+1:end) = nx*(ny-1)+1:nx*ny;

elseif P == 2 %Eriksson-Johnson problem 
    dir = zeros(1,ny+2*nx-2);
    dir(1:nx) = 1:nx;
    for i = 1:ny-2
        dir(nx+2*i-1) = i*nx+1;
        dir(nx+i) = (i+1)*nx;
    end
    dir(nx+ny-1:end) = nx*(ny-1)+1:nx*ny;

elseif P == 3 %Hochbruck-Osterman equation
    dir = zeros(1,2*ny+2*nx-4);
    dir(1:nx) = 1:nx;
    for i = 1:ny-2
       dir(nx+2*i-1) = i*nx+1;
       dir(nx+2*i) = (i+1)*nx;
    end
    dir(nx+2*(ny-2)+1:end) = nx*(ny-1)+1:nx*ny;
    
elseif P == 4 %Allen-Cahn equation
    dir=[];
end

end