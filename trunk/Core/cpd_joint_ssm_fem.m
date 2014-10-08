function [TY1, TY2, R1, R2, t1, t2, s1, s2, b, fem, U1, U2]= cpd_joint_ssm_fem(X1, X2, SSM, w, errtol, maxiters, R1, R2, t1, t2, s1, s2, nMods, b, sigma21, sigma22, mu, beta, E, nu)
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

D = size(X1,2);
N1 = size(X1,1);
N2 = size(X2,1);

M = size(SSM.mean,1);

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
if (nargin < 7 || isempty(R1))
    R1 = eye(D);
end
if (nargin < 8 || isempty(R2))
    R2 = eye(D);
end
if (nargin < 9 || isempty(t1))
    t1 = zeros(D, 1);
end
if (nargin < 10 || isempty(t2))
    t2 = zeros(D, 1);
end
if (nargin < 11 || isempty(s1))
    s1 = 1;
end
if (nargin < 12 || isempty(s2))
    s2 = 1;
end
if (nargin < 13 || isempty(nMods))
    nMods = 5;
end
if (nargin < 14 || isempty(b))
    b = zeros(nMods,1);
end

% Initialize the registration
Psi = SSM.mods(:,1:nMods);

z_vec = reshape(SSM.mean',[],1);

y_vec = z_vec+Psi*b;

v_vec1 = zeros(D*M,1);
v_vec2 = zeros(D*M,1);

Lambda = SSM.latent(1:nMods);
Lambda = ones(nMods,1)./Lambda;
Lambda = spdiags(Lambda,0,nMods,nMods);

Y = reshape(y_vec,D,[])';
TY1 = reshape(y_vec+v_vec1,D,[])';
TY1 = bsxfun(@plus, TY1*(R1')*s1, t1');

TY2 = reshape(y_vec+v_vec2,D,[])';
TY2 = bsxfun(@plus, TY2*(R2')*s2, t2');

if (nargin < 15 || isempty(sigma21))
    % estimate initial variance
    XX = reshape(X1, [1, N1, D]);
    YY = reshape(TY1, [M, 1, D]);
    XX = repmat(XX, [M, 1, 1]);
    YY = repmat(YY, [1, N1, 1]);
    diff = XX-YY;
    diff = diff.*diff;
    err2 = sum(diff(:));
    sigma21 = 1/(D*N1*M)*err2;
    clear diff;
end
if (nargin < 16 || isempty(sigma22))
    XX = reshape(X2, [1, N2, D]);
    YY = reshape(TY2, [M, 1, D]);
    XX = repmat(XX, [M, 1, 1]);
    YY = repmat(YY, [1, N2, 1]);
    diff = XX-YY;
    diff = diff.*diff;
    err2 = sum(diff(:));
    sigma22 = 1/(D*N2*M)*err2;
    clear diff;
end

% FEM stuff
D_material = fem_material_linear(E, nu); 
[nodes, ~, elems] = tetgen_mex(Y', (SSM.faces)',[],'');
Phi = [speye(M),zeros(M, size(nodes,2)-M)];
Phi_tilde = kron(Phi, eye(D));
fem = fem_model(nodes', elems');
setAbortJacobian(fem, EPSILON);      % abort on inverted elements

iters = 0;
err1 = errtol+1;
err2 = errtol+1;

C = ones(D,1);  % temp, so we don't keep recreating, for rigid transform

while ((iters < maxiters) && (err1 > errtol || err2 > errtol ) && (sigma21 > 0.1 || sigma22 > 0.1))
    
    % E-step
    [P1, P11, Pt11, Np1] = cpd_P(X1, TY1, sigma21, w);
    
    [P2, P12, Pt12, Np2] = cpd_P(X2, TY2, sigma22, w);
    
    % Here, TY is just Y+V
    TY1 = y_vec+v_vec1;
    TY1 = reshape(TY1,D,[])';
    
    TY2 = y_vec+v_vec2;
    TY2 = reshape(TY2,D,[])';
    
    % M-step
    % Rigid align
    mux1 = sum(P1*X1,1)/Np1;
    muy1 = sum(P1'*TY1,1)/Np1;
    
    mux2 = sum(P2*X2,1)/Np2;
    muy2 = sum(P2'*TY2,1)/Np2;

    % recompute because P changed, causing mux/y to change
    XX = bsxfun(@minus, X1, mux1);
    YY = bsxfun(@minus, TY1, muy1);
    A = XX'*P1'*YY;    
    [U, ~, V] = svd(A);
    C(D) = det(U*V');
    R1 = U*diag(C)*V';
    %s1 = trace(A'*R1)/trace(YY'*diag(P11)*YY)
    s1 = 1.0;
    t1 = mux1' - s1*R1*(muy1');
    
    XX = bsxfun(@minus, X2, mux2);
    YY = bsxfun(@minus, TY2, muy2);
    A = XX'*P2'*YY;    
    [U, ~, V] = svd(A);
    C(D) = det(U*V');
    R2 = U*diag(C)*V';
    %s2 = trace(A'*R2)/trace(YY'*diag(P12)*YY)
    s2 = 1.0;
    t2 = mux2' - s2*R2*(muy2');

    % Shape update
    dP1 = spdiags(kron(P11,ones(D,1)),0,D*M,D*M);
    dP2 = spdiags(kron(P12,ones(D,1)),0,D*M,D*M);
    LHS = s1*s1*sigma22*sigma22*Psi'*dP1*Psi + s2*s2*sigma21*sigma21*Psi'*dP2*Psi + mu*sigma21*sigma22*Lambda;
    
    RHS1 = s1*(P1*X1)*R1;
    RHS1 = sigma22*Psi'*(reshape(RHS1',[],1) - dP1*(s1*s1*(z_vec+v_vec1)+s1*repmat(R1'*t1,M,1)));
    
    RHS2 = s2*(P2*X2)*R2;
    RHS2 = sigma21*Psi'*(reshape(RHS2',[],1) - dP2*(s2*s2*(z_vec+v_vec2)+s2*repmat(R2'*t2,M,1)));
    
    RHS = RHS1 + RHS2;
    
    b = LHS\RHS;
    
    % Change model, update stiffness
    y_vec = z_vec+Psi*b;
    Y = reshape(y_vec,D,[])';
    
    % update nodes in existing model
    % model is independent of rotation/scale... in SSM coordinates
    
        updateNodes(fem,Y,1:M);
        
        % compute stiffness
        [K, minJ] = getStiffnessMatrix(fem, D_material);
        if (isempty(K))
            fprintf('Inverted element detected (%g), regenerating mesh... ',...
                minJ);
            tic;
            [nodes, ~, elems] = tetgen_mex(Y', (SSM.faces)', 2, '');
            rtime = toc;
            fprintf('%g s\n', rtime);
            set(fem, nodes', elems');
            K = getStiffnessMatrix(fem, D_material);
            Phi = [speye(M),zeros(M, size(nodes,2)-M)];
            Phi_tilde = kron(Phi, eye(D));
        else
            %fprintf('Mesh can be re-used, minJ=%f\n',minJ);
        end
    
        % FEM step
        LHS = s1*s1*(Phi_tilde')*dP1*Phi_tilde + beta*sigma21*K;
        RHS = -s1*(P1*X1)*R1;
        RHS = -Phi_tilde'*(reshape(RHS',[],1)+dP1*(s1*s1*y_vec+s1*repmat(R1'*t1,M,1)));
        u_vec1 = LHS\RHS;
        v_vec1 = Phi_tilde*u_vec1;
        
        LHS = s2*s2*(Phi_tilde')*dP2*Phi_tilde + beta*sigma22*K;
        RHS = -s2*(P2*X2)*R2;
        RHS = -Phi_tilde'*(reshape(RHS',[],1)+dP2*(s2*s2*y_vec+s2*repmat(R2'*t2,M,1)));
        u_vec2 = LHS\RHS;
        v_vec2 = Phi_tilde*u_vec2;
    
    % Update GMM centroids
    TY1 = y_vec+v_vec1;
    TY1 = reshape(TY1,D,[])';
    TY1 = bsxfun(@plus, TY1*(R1')*s1, t1');
    
    TY2 = y_vec+v_vec2;
    TY2 = reshape(TY2,D,[])';
    TY2 = bsxfun(@plus, TY2*(R2')*s2, t2');
    
    % Update sigma
    PX1 = P1*X1;
    xPx = (Pt11')*sum(X1.*X1,2);
    yPy = (P11')*sum(TY1.*TY1,2);
    trPXTY = sum(TY1(:).*PX1(:));
    
    err1 = sigma21;
    sigma21 = (xPx-2*trPXTY+yPy)/(Np1*D);
    
    PX2 = P2*X2;
    xPx = (Pt12')*sum(X2.*X2,2);
    yPy = (P12')*sum(TY2.*TY2,2);
    trPXTY = sum(TY2(:).*PX2(:));
    
    err2 = sigma22;
    sigma22 = (xPx-2*trPXTY+yPy)/(Np2*D);
    
    err1 = abs((err1-sigma21)/err1);	% percent change
    err2 = abs((err2-sigma22)/err2);	% percent change
    
    disp(['Error fraction:', num2str(err1), ' ', num2str(err2)]);
    
    TYstiff = y_vec;
    TYstiff = reshape(TYstiff,D,[])';
    
    [az, el] = view;
    clf;
    
    subplot(1,2,1); axis([-5,5,-5,5,-5,5]); plot3(X1(:,1), X1(:,2), X1(:,3),'.r','MarkerSize',10); hold on;
    
    % without deformation
    
    TYstiff1 = bsxfun(@plus, TYstiff*(R1')*s1, t1');
    plot3(TYstiff1(:,1), TYstiff1(:,2), TYstiff1(:,3),'og','MarkerSize', 5, 'MarkerFaceColor', 'g');
    plot3(TY1(:,1), TY1(:,2), TY1(:,3),'ob','MarkerSize', 5, 'MarkerFaceColor', 'b');
    legend({'Target','SSM only','SSM-FEM'},'location','North');
    hold off;
    
    subplot(1,2,2); axis([-5,5,-5,5,-5,5]); plot3(X2(:,1), X2(:,2), X2(:,3),'.r','MarkerSize',10); hold on;
    
    % without deformation
    TYstiff2 = bsxfun(@plus, TYstiff*(R2')*s2, t2');
    plot3(TYstiff2(:,1), TYstiff2(:,2), TYstiff2(:,3),'og','MarkerSize', 5, 'MarkerFaceColor', 'g');
    plot3(TY2(:,1), TY2(:,2), TY2(:,3),'ob','MarkerSize', 5, 'MarkerFaceColor', 'b');
    legend({'Target','SSM only','SSM-FEM'},'location','North');

    hold off;
    
    set(gcf,'color','w');
    drawnow;
    
    %saveas(gcf,num2str(iters),'png');
    export_fig( num2str(iters), '-png');
        
    iters = iters + 1;
end

U1 = reshape(u_vec1,D,[])';
U2 = reshape(u_vec2,D,[])';

end