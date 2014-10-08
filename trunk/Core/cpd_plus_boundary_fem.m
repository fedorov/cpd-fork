function [TY, V, sigma2, fem, U, Phi] = cpd_plus_boundary_fem(X, Y, F, lambda, beta2, w, errtol, maxiters, sigma2, E, nu )
% Does deformable CPD, then uses the new boundary and a FEM to interpolate
% interior points.

D = size(X,2);
N = size(X,1);
M = size(Y,1);

% regular CPD
[TY, P, sigma2] = cpd_coherent(X, Y, lambda, beta2, w, errtol, maxiters, sigma2);

% create original FEM
[nodes, ~, elems] = tetgen_mex(Y', F',[],'');
fem = fem_model(nodes', elems');
Nfem = size(nodes,2);

% FEM stuff
D_material = fem_material_linear(E, nu); 

% get stiffness matrix
[K, ~] = getStiffnessMatrix(fem, D_material);

% get interpolation matrix
Phi = [speye(M),zeros(M, size(nodes,2)-M)];

% set boundary conditions
K(1:(3*M), 1:(3*M)) = speye(3*M, 3*M);
b = zeros(3*Nfem, 1);
b(1:(3*M)) = TY(:)-Y(:);

% solve for displacements
u_vec = K\b;

U = reshape(u_vec,D,[])';
V = Phi*U;

end