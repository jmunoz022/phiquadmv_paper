function [dirx,diry] = BC_1D(nx,ny,P)
%
% BC_1D - Function that returns Dirichlet nodes in the 1D meshes 
%      for each problem
%                   
% INPUT: 
%   nx - number of nodes in the x direction
%   ny - number of nodes in the x direction
%   P  - Problem number:
%        - P=1: Heat equation
%        - P=2: Eriksson-Johnson problem 
%        - P=4: %Hochbruck-Osterman equation
%        - P=3: Allen-Cahn equation 
%
% OUTPUT: 
%   dirx - vector containing the Dirichlet nodes in the x direction
%   diry - vector containing the Dirichlet nodes in the x direction

if P == 1 %Heat Equation
    dirx = [1 nx];
    diry = [1 ny];
    
elseif P == 2 %Eriksson-Johnson problem
    dirx = nx;
    diry = [1 ny];

elseif P == 3 %Hochbruck-Osterman equation
    dirx = [1 nx];
    diry = [1 ny];
    
elseif P == 4 %Allen-Cahn equation
    dirx = [];
    diry = [];
    
end

end