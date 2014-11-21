function [TY, R, t, s, fem, u, sigma2, P] = cpd_rigid_fem(X, fX, Y, fY, w, errtol, maxiters, R, t, s, sigma2, beta, E, nu)
%cpd_biomech performs biomechanically constrained point cloud registration
%without an SSM
%   [TY, sigma2] = cpd_biomech(X, fX, Y, fY, w, errtol, maxiters, R, t, s, sigma2, beta, E, nu)
%
%   Inputs:
%   X    N x D matrix of input points
%   fX   matrix of input faces corresponding to input points (optional)
%   Y    M x D matrix of Gaussian centres
%   fY   matrix faces corresponding to Gaussian centres
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
%   sigma2  Initial estimate of squared variance (optional, default
%        estimated from X, Y)
%   beta    Weight to incorporate FEM
%   E,nu    Young's modulus and Poisson's ratio
%
%   Outputs:
%   TY   Transformed points%   (R, t, s)  Final transformation
%   P    Alignment probabilities
%   sigma2  Estimated probability variance

EPSILON = 1e-16;	% default small value

D = size(X,2);
N = size(X,1);
M = size(Y,1);

% allow variable input args
if (nargin < 5 || isempty(w))
    w = 0;
end
if (nargin < 6 || isempty(errtol))
    errtol = 1e-5;
end
if (nargin < 7 || isempty(maxiters))
    maxiters = 100;
end
if (nargin < 8 || isempty(R))
    R = eye(D);
end
if (nargin < 9 || isempty(t))
    t = zeros(D, 1);
end
if (nargin < 10 || isempty(s))
    s = 1;
end

% Initialize the registration
TY = bsxfun(@plus, Y*(R')*s, t');
y_vec = reshape(TY',[],1);
v_vec = zeros(D*M,1);

if (nargin < 11 || isempty(sigma2))
    % estimate initial variance
    XX = reshape(X, [1, N, D]);
    YY = reshape(TY, [M, 1, D]);
    XX = repmat(XX, [M, 1, 1]);
    YY = repmat(YY, [1, N, 1]);
    diff = XX-YY;
    diff = diff.*diff;
    err2 = sum(diff(:));
    sigma2 = 1/(D*N*M)*err2;
    clear diff;
end

% FEM stuff
D_material = fem_material_linear.getElasticity.getElasticity(E, nu); 
[nodes, ~, elems] = tetgen_mex(Y', (fY)',[],'');
Phi = [speye(M),zeros(M, size(nodes,2)-M)];
Phi_tilde = kron(Phi, eye(D));
fem = fem_model(nodes', elems');
[K, minJ] = getStiffnessMatrix(fem, D_material);
setAbortJacobian(fem, EPSILON);      % abort on inverted elements

iters = 0;
err = errtol+1;

C = ones(D,1);  % temp, so we don't keep recreating, for rigid transform    && (err > errtol)
while ((iters < maxiters) )
    
    % E-step
    [P, P1, Pt1, Np] = cpd_P(X, TY, sigma2, w);
    
    % Here, TY is just Y+V
    TY = y_vec+v_vec;
    TY = reshape(TY,D,[])';
    
    % M-step
    % Rigid align
    mux = sum(P*X,1)/Np;
    muy = sum(P'*TY,1)/Np;

    % recompute because P changed, causing mux/y to change
    XX = bsxfun(@minus, X, mux);
    YY = bsxfun(@minus, TY, muy);
    A = XX'*P'*YY;    
    [U, ~, V] = svd(A);
    C(D) = det(U*V');
    R = U*diag(C)*V';
    s=1;
    %s = trace(A'*R)/trace(YY'*diag(P1)*YY);
    t = mux' - s*R*(muy');

    % FEM step
    dP1 = spdiags(kron(P1,ones(D,1)),0,D*M,D*M);
    LHS = s*s*(Phi_tilde')*dP1*Phi_tilde + beta*sigma2*K;
    RHS = -s*(P*X)*R;
    RHS = -Phi_tilde'*(reshape(RHS',[],1)+dP1*(s*s*y_vec+s*repmat(R'*t,M,1)));
    u_vec = LHS\RHS;
    v_vec = Phi_tilde*u_vec;
    
    % Update GMM centroids
    TY = y_vec+v_vec;
    TY = reshape(TY,D,[])';
    TY = bsxfun(@plus, TY*(R')*s, t');
    
    % Update sigma
    PX = P*X;
    xPx = (Pt1')*sum(X.*X,2);
    yPy = (P1')*sum(TY.*TY,2);
    trPXTY = sum(TY(:).*PX(:));
    
    err = sigma2;
    sigma2 = (xPx-2*trPXTY+yPy)/(Np*D);
    err = (err-sigma2)/err;		% percent change
    
%     clf;
%     patch('Vertices',X,'Faces',fX,'FaceColor','r','FaceAlpha',.5);
%     hold on;
%     patch('Vertices',TY,'Faces',fY,'FaceColor','b','FaceAlpha',.5);
%     view(gca,-88,22);
%     drawnow;
%     
%     saveas(gcf,num2str(iters),'png');
    
    u = reshape(u_vec,D,[])';
    
    iters = iters + 1;
end

V = reshape(v_vec,D,[])';

end

