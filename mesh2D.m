function [xsol,ysol,nelx,nely] = mesh2D(r,x1,x2,y1,y2,P)
%
% MESH2D - Function that returns the discretization in space and the number of elements
%          in each space direction (2^r)
%
% INPUT: 
%   r             - exponent of 2 
%   [x1 x2 y1 y2] - limits of the rectangular domain in space
%   P             - problem number:
%                   - P=1: Heat equation
%                   - P=2: Eriksson-Johnson problem 
%                   - P=3: Hochburck-Osterman equation 
%                   - P=4: Allen-Cahn equation  
%
% OUTPUT: 
%   xsol - spatial grid in the x direction         
%   ysol - spatial grid in the y direction
%   nelx - number of elements in the x direction
%   nely - numer of elements in the y direction
%
if P == 1 %Heat equation: Regular mesh
    nelx = 2^r;
    nely = 2^r;
    xsol = linspace(x1,x2,nelx+1);
    ysol = linspace(y1,y2,nely+1);

elseif P == 2 %Eriksson-Johnson problem: Shiskin mesh
    nelx = 2^r;
    nely = 2^r;
    xsol0 = linspace(x1,x2-(x2-x1)/16,nelx/2+1);
    xsol1 = linspace(x2-(x2-x1)/16,x2,nelx/2+1);
    xsol = union(xsol0,xsol1);
    ysol = linspace(y1,y2,nely+1);

elseif P == 3 %Hochburck-Osterman equation: Regular mesh 
    nelx = 2^r;
    nely = 2^r;
    xsol = linspace(x1,x2,nelx+1);
    ysol = linspace(y1,y2,nely+1);
    
elseif P == 4 %Allen-Cahn equation: Regular mesh excluding x=y=0
    nelx = 2^r+1;
    nely = 2^r+1;
    xsol = linspace(x1,x2,nelx+1);
    ysol = linspace(y1,y2,nely+1);
end

end

