function [Ax,Ay,A,U0,M] = compute_matrices(P,r,steps)
%
% COMPUTE_MATRICES - Function that returns the full 2D sparse matrix A,
%                    the 1D matrices Ax and Ay, and the initial 
%                    condition U0
%
% INPUT: 
%   P - problem number:
%           - P=1 for the Heat equation
%           - P=2 for the Eriksson-Johnson probem
%   r - number of elements in each space dimension: 2^r
%   steps - number of time steps: T/steps
%
% OUTPUT: 
%   Ax,Ay - the 1D sparse matrices
%   A - the full 2D sparse matrix
%   U0 - the initial condition

%Import the data of the problem
[eps,beta,u0,x1,x2,y1,y2,T] = data(P);

%Meshes and Parameters
[xsol,ysol,nelx,nely] = mesh2D(r,x1,x2,y1,y2,P);
[nx,ny,nel,nnode,coord,nodes] = parameters(xsol,ysol);
tau = T/steps; %time step size

%1D matrices
nquad = 2; %quadrature points
[dirx,diry] = BC_1D(nx,ny,P);
Ax = -tau*sparse(FEM_matrices_1D(eps,beta(1),dirx,xsol,nelx,nquad));
Ay = -tau*sparse(FEM_matrices_1D(eps,beta(2),diry,ysol,nely,nquad));

%Full sparse matrix
dir = BC(nx,ny,P);
[A,M] = FEM_matrices_sparse(eps,beta,dir,nel,coord,nodes,nquad);
A = -tau*A;
%Akron=kron(eye(size(Ay)),Ax)+kron(Ay,eye(size(Ax))); norm(full(A-Akron))

%Initial condition vector
[xmat,ymat] = meshgrid(xsol,ysol);
U0 = u0(xmat,ymat);
U0 = reshape(U0',[],1);
U0(dir) = [];

end