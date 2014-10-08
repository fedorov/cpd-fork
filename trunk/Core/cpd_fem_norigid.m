function [TY, V, sigma2, P, fem, U, Phi] = cpd_fem_norigid(X, Y, F, w, errtol, maxiters, sigma2, beta, E, nu, ...
    fem, Phi, FV)
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
if (nargin < 4 || isempty(w))
    w = 0;
end
if (nargin < 5 || isempty(errtol))
    errtol = 1e-5;
end
if (nargin < 6 || isempty(maxiters))
    maxiters = 100;
end

% Initialize the registration
y_vec = reshape(Y',[],1);
v_vec = zeros(D*M,1);

Y = reshape(y_vec,D,[])';
TY = reshape(y_vec+v_vec,D,[])';

if (nargin < 7 || isempty(sigma2))
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

if (nargin < 11 || isempty(fem))
    [nodes, ~, elems] = tetgen_mex(Y', F',[],'');
    fem = fem_model(nodes', elems');
    setAbortJacobian(fem, EPSILON);      % abort on inverted elements
end

if (nargin < 12 || isempty(Phi))
    Phi = [speye(M),zeros(M, size(nodes,2)-M)];
end

if (nargin < 13 || isempty(FV))
    FV = [];
end

% FEM stuff
D_material = fem_material_linear(E, nu); 
% compute stiffness
[K, minJ] = getStiffnessMatrix(fem, D_material);
if (minJ < 0) 
    error('Stiffness matrix has inverted elements');
end
Phi_tilde = kron(Phi, eye(D));

iters = 0;
err = errtol+1;

C = ones(D,1);  % temp, so we don't keep recreating, for rigid transform    
while ((iters < maxiters) && (err > errtol))
    
    % E-step
    [P, P1, Pt1, Np] = cpd_P(X, TY, sigma2, w);
    P1(P1<1e-10) = 1e-10;
    
    % Here, TY is just Y+V
    TY = y_vec+v_vec;
    TY = reshape(TY,D,[])';
    
    % M-step
    % FEM step
    dP1 = spdiags(kron(P1,ones(D,1)),0,D*M,D*M);
    LHS = (Phi_tilde')*dP1*Phi_tilde + beta*sigma2*K;
    RHS = -(P*X);
    RHS = -Phi_tilde'*(reshape(RHS',[],1)+dP1*(y_vec));
    u_vec = LHS\RHS;
    v_vec = Phi_tilde*u_vec;
    
    % Update GMM centroids
    TY = y_vec+v_vec;
    TY = reshape(TY,D,[])';
    
    % Update sigma
    PX = P*X;
    xPx = (Pt1')*sum(X.*X,2);
    yPy = (P1')*sum(TY.*TY,2);
    trPXTY = sum(TY(:).*PX(:));
    
    err = sigma2;
    sigma2 = (xPx-2*trPXTY+yPy)/(Np*D);
    err = abs((err-sigma2)/err);	% percent change
    disp(['Error fraction:', num2str(err)]);
    
        [az, el] = view;
        clf;
        if (isempty(FV))
            plot3(X(:,1), X(:,2), X(:,3),'.r','MarkerSize',10);
        else 
            patch(FV,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
        end
        hold on;
        % axis([-3,3,-3,3,-3,3]);
        plot3(TY(:,1), TY(:,2), TY(:,3),'ob','MarkerSize', 5, 'MarkerFaceColor', 'b');
        % axis([-3,3,-3,3,-3,3]);
        legend({'target','fem'},'location','NEo')
        view(az, el);
        drawnow;
        
    iters = iters + 1;
end

V = reshape(v_vec,D,[])';
U = reshape(u_vec,D,[])';

end