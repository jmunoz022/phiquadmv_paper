function [A,M] = FEM_matrices_1D(eps,beta,dir,sol,nel,nquad)
%
% FEM_MATRICES_1D - Function that returns the 1D matrix
%                   coming from the spatial discretization
%
% INPUT: 
%   eps   - diffusion coefficient 
%   beta  - advection vector
%   dir   - Dirichlet nodes
%   sol   - grid in space
%   nel   - total number of elements
%   quad  - number of quadrature points 
%
% OUTPUT: 
%   A - the 1D matrix
%

%Total number of nodes
nnode = size(sol,2);

%Numbering of nodes
nodes = zeros(nel,2);
for i = 1:nel
    nodes(i,1) = i;
    nodes(i,2) = i+1;
end

%Initialization of the matrices
kk = zeros(nnode,nnode); 
mm = zeros(nnode,nnode); 
gg = zeros(nnode,nnode); 

%Loop through elements
for iel = 1:nel
    
    %Get node and coordinates
    nd = nodes(iel,:);
    xcoord = sol(nd);
    
    %Compute local matrices
    m = local_mass_1D(xcoord,nquad);
    k = local_stiffness_1D(xcoord,eps,nquad);
    g = local_advection_1D(xcoord,beta,nquad);
    
    %Assembly
    mm(nd,nd) = mm(nd,nd)+m;
    kk(nd,nd) = kk(nd,nd)+k;
    gg(nd,nd) = gg(nd,nd)+g;
end
%Transpose matrices 
M = sparse(mm');
K = sparse(kk');
G = sparse(gg');

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