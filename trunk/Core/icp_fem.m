function [TY, V, fem, U, Phi] = icp_fem(X, Y, F, errtol, maxiters, R, t, s, beta, E, nu, fem, Phi, FV)
%cpd_biomech performs biomechanically constrained point cloud registration
%   [TY, sigma2] = cpd_biomech(X, Y, w, errtol, maxiters, R, t, s, sigma2, beta )
%
%   Inputs:
%   X    N x D matrix of input points
%   Y    M x D matrix of Gaussian centres
%   w    Weight to account for noise and outliers (optional, defaults to 0)
%   errtol Error tolerance for convergence.  Algorithm will terminate if
%        the objective function does not change by more than this amount
%        (optional, defaults to 1e-10)
%   maxiters Maximum number of iterations.  Algorithm will terminate if
%        this many iterations have been performed (optional, defaults to
%        100)
%   R    Initial D x D affine matrix (optional, defaults to identity)
%   t    Initial D x 1 translation vector (optional, defaults to zeroes)
%   s    single scale value
%   beta    Weight to incorporate FEM
%   sigma2  Initial estimate of squared variance (optional, default
%        estimated from X, Y)
%
%   Outputs:
%   TY   Transformed points
%   (B, t)  Final transformation
%   P    Alignment probabilities
%   sigma2  Estimated probability variance

EPSILON = 1e-16;	% default small value

D = size(X,2);
N = size(X,1);
M = size(Y,1);

% allow variable input args
if (nargin < 4 || isempty(errtol))
    errtol = 1e-5;
end
if (nargin < 5 || isempty(maxiters))
    maxiters = 100;
end
if (nargin < 6 || isempty(R))
    R = eye(D);
end
if (nargin < 7 || isempty(t))
    t = zeros(D, 1);
end
if (nargin < 8 || isempty(s))
    s = 1;
end

% Initialize the registration
y_vec = reshape(Y',[],1);
v_vec = zeros(D*M,1);

Y = reshape(y_vec,D,[])';
TY = reshape(y_vec+v_vec,D,[])';
TY = bsxfun(@plus, TY*(R')*s, t');

if (nargin < 9 || isempty(fem))
    [nodes, ~, elems] = tetgen_mex(Y', F',[],'');
    fem = fem_model(nodes', elems');
    setAbortJacobian(fem, EPSILON);      % abort on inverted elements
end

if (nargin < 10 || isempty(Phi))
    Phi = [speye(M),zeros(M, size(nodes,2)-M)];
end

if (nargin < 11 || isempty(FV))
    FV = [];
end

% FEM stuff
D_material = fem_material_linear.getElasticity(E, nu);

% compute stiffness
[K, minJ] = getStiffnessMatrix(fem, D_material);
if (minJ < 0)
    error('Stiffness matrix has inverted elements');
end

iters = 0;
err = errtol+1;

u_vec = [];

Phi_tilde = kron(Phi, eye(D));

% Xtree = kdtree_build( X );

while ((iters < maxiters) && (err > errtol))
    
    % Find correspondences
    tree = kdtree_build( TY );
    idxs = kdtree_nearest_neighbor(tree, X);
    kdtree_delete(tree);
    clear tree;
    
    % Find surface forces
    Fs = zeros(size(TY));
    Fs(idxs,:) = X - TY(idxs,:);
    
    % add forces from the other side
    % idxs = kdtree_nearest_neighbor(Xtree, TY);
    % Fs = 0.5*(Fs + X(idxs,:)-TY);
    
    Fs_vec = reshape(Fs',[],1);
    
    % for i=1:size(Fs,1)
    %    plot3([TY(i,1); TY(i,1)+Fs(i,1)], [TY(i,2); TY(i,2)+Fs(i,2)], [TY(i,3); TY(i,3)+Fs(i,3)], 'g', 'LineWidth',1);
    % end
    
    % Solve linear system
    s1 = warning('error', 'MATLAB:singularMatrix');
    warning('error', 'MATLAB:singularMatrix');
    if ( ~isempty(u_vec) )
        u_vec_old = u_vec;
    else
        u_vec_old = zeros(size(K,1),1);
    end
    
    LHS = speye(size(K,1)) + beta*K;
    RHS = u_vec_old + beta*(Phi_tilde')*Fs_vec;
    
    try
        % Regular processing part
        u_vec = LHS\RHS;
    catch
        % Exception-handling part
        fprintf('Can''t solve linear system (reason: %s)\n', lasterr);
        u_vec = u_vec_old;
        iters = maxiters;
    end
    warning(s1);
    
    v_vec = Phi_tilde*u_vec;
    
    TY = y_vec+v_vec;
    TY = reshape(TY,D,[])';
    
    [az, el] = view;
    clf;
    if (isempty(FV))
        plot3(X(:,1), X(:,2), X(:,3),'.r','MarkerSize',10);
    else
        patch(FV,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
    end
    
    hold on;
    patch('Vertices', TY, 'Faces', F, 'FaceColor','blue','FaceAlpha',0.2, 'EdgeColor', 'blue', 'EdgeAlpha',0.2);
    legend({'target','fem'},'location','NEo')
    view(az, el);
    drawnow;
    % saveas(gcf,num2str(iters),'png');
    
    iters = iters + 1;
end

V = reshape(v_vec,D,[])';
U = reshape(u_vec,D,[])';

% kdtree_delete(Xtree);
% clear Xtree

end