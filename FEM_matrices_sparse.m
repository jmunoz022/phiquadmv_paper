function [A,M] = FEM_matrices_sparse(eps,beta,dir,nel,coord,nodes,nquad)
%
% FEM_MATRICES_SPARSE - Function that returns the full 2D sparse matrix
%                       coming from the spatial discretization
%
% INPUT: 
%   eps   - diffusion coefficient 
%   beta  - advection vector
%   dir   - Dirichlet nodes
%   nel   - total number of elements
%   coord - vector containing the physical coordinates
%   nodes - vector numbering the nodes
%   quad  - number of quadrature points 
%
% OUTPUT: 
%   A - the 2D sparse matrix M^-1(G+K)
%   M - the 2D sparse mass matrix 
%

%Initialization of indexes and values
I = zeros(16*nel,1);
J = zeros(16*nel,1);
Vk = zeros(16*nel,1);
Vm = zeros(16*nel,1);
Vg = zeros(16*nel,1);

%Loop through elements 
for iel = 1:nel
    
    %Get node and coordinates
    nd = nodes(iel,:);
    xcoord = coord(nd,1);
    ycoord = coord(nd,2);
    
    %Local matrices
    k = local_stiffness(xcoord,ycoord,eps,nquad);
    g = local_advection(xcoord,ycoord,beta,nquad);
    m = local_mass(xcoord,ycoord,nquad);
    
    
    %Indexes
    I((iel-1)*16+1:iel*16) = repmat(nd',4,1);
    Jaux = repmat(nd,4,1);
    J((iel-1)*16+1:iel*16) = Jaux(:);
    
    %Values
    Vk((iel-1)*16+1:iel*16) = k(:);
    Vm((iel-1)*16+1:iel*16) = m(:);
    Vg((iel-1)*16+1:iel*16) = g(:);
    
end
%Assemble the matrices
kk = sparse(I,J,Vk);
mm = sparse(I,J,Vm);
gg = sparse(I,J,Vg);

%Transpose the matrices
M = mm';
K = kk';
G = gg';

%Eliminate Dirichlet degrees of freedom
M(dir,:) = [];
M(:,dir) = [];
K(dir,:) = [];
K(:,dir) = [];
G(dir,:) = [];
G(:,dir) = [];

%Calculate matrix A
A = M\(G+K);

end