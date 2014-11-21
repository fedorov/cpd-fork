function [ TY, P, sigma2 ] = tps_rpm( X, Y, w, errtol, maxiters, lambda1, lambda2, annealing, sigma2)
%CPD_COHERENT performs the coherent CPD algorithm
%   [ TY, P, sigma2 ] = cpd_coherent( X, Y, lambda, beta, w, errtol, maxiters, sigma2)
%
%   Inputs:
%   X    N x D matrix of input points
%   Y    M x D matrix of Gaussian centres
%   lambda > 0, regularization weight
%   beta2  > 0, factor controlling amount of coherence
%   w    Weight to account for noise and outliers (optional, defaults to 0)
%   errtol Error tolerance for convergence.  Algorithm will terminate if
%        the objective function does not change by more than this amount
%        (optional, defaults to 1e-10)
%   maxiters Maximum number of iterations.  Algorithm will terminate if
%        this many iterations have been performed (optional, defaults to
%        100)
%   sigma2  Initial estimate of squared variance (optional, default
%        estimated from X, Y)
%   
%   Ouputs:
%   TY   Transformed Gaussian points.  The transform can then be recovered
%        by subtracting this from the original points
%   P    Alignment probabilities
%   sigma2 Estimated probability variance

    D = size(X,2);
    N = size(X,1);
    M = size(Y,1);
    
    % allow variable input args
    if (nargin < 5 || isempty(w))
        w = 0;
    end
    if (nargin < 6 || isempty(errtol))
        errtol = 1e-10;
    end
    if (nargin < 7 || isempty(maxiters))
        maxiters = 100;
    end
    
    % initial transformed input
    TY = Y;
    
    if (nargin < 8 || isempty(sigma2))
        % estimate initial variance
        % Restored to remove additional function call
        XX = reshape(X, [1, N, D]);
        YY = reshape(TY, [M, 1, D]);
        XX = repmat(XX, [M, 1, 1]);
        YY = repmat(YY, [1, N, 1]);
        diff = XX-YY;
        diff = diff.*diff;
        err2 = sum(diff(:));
        sigma2 = 1/(D*N*M)*err2;
    end
    
    % G matrix
    XX = reshape(Y, [1, M, D]);
    YY = reshape(Y, [M, 1, D]);
    XX = repmat(XX, [M, 1, 1]);
    YY = repmat(YY, [1, M, 1]);
    diff = XX-YY;
    G = sqrt(sum(diff.*diff,3));
    clear diff XX YY;
    
    % Transform parameters
    % W = zeros(M,D);
    
    % initialize loop
    iters = 0;
    err = errtol+1;  % initialize so we enter loop
    % q = -Inf;   % force first q to be > errtol away from next
    W = [];
    perT_maxit = 5;
    while ((iters < maxiters) && (err > errtol))
        
        for i=1:perT_maxit
            % E-step
            [P, P1, Pt1, Np] = cpd_P(X, TY, sigma2, w);
            
            T0 = max(TY(:,1))^2;
            moutlier = 1/sqrt(T0)*exp(-1);
            m_outliers_row = ones(1,N) * moutlier; 
            m_outliers_col = ones(M,1) * moutlier; 
            
            m = cMIX_calc_m (vx, y, T, m_method, m_outliers_row, m_outliers_col,it_total,k_sigma);
            
            % Given v=tranformed(x), update m:
            y_tmp = zeros (xmax, ymax);
            for it_dim=1:dim
                y_tmp = y_tmp + (vx(:,it_dim) * ones(1,ymax) - ones(xmax,1) * y(:,it_dim)').^2;
            end;
            
            m_tmp = exp (-y_tmp/T);
            m_tmp = m_tmp + randn(xmax, ymax) * (1/xmax) * 0.001;
            
            % double normalization, but keep outlier entries constant.
            moutlier       = 1/xmax * 0.1;
            m_outliers_row = ones (1,ymax) * moutlier;
            m_outliers_col = ones (xmax,1) * moutlier;
            
            [m, junk1, junk2] = cMIX_normalize_m (m_tmp, m_outliers_col, m_outliers_row);
            
            % M-step
            P1(P1<1e-10) = 1e-10;
            P1i = 1./P1;
            PX = P*X;
            
            v = [ones(M,1), TY];
            
            [q,r] = qr(v);
            q1 = q(:,1:D+1);
            q2 = q(:,D+2:M);
            R = r(1:D+1,1:D+1);
            
            x = [ones(M,1), TY-PX];
            
            LHS = q2'*G*q2 + lambda1*eye(M-(D+1));
            RHS = q2'*x;
            gamma = LHS\RHS;
            W = q2*gamma;
            
            LHS = R'*R + lambda2 * eye(length(R),length(R));
            RHS = R'*q1'*(x-G*q2*gamma) - R'*R;
            A =  LHS\RHS;
            d = A + eye(D+1,D+1)
            trace((d-eye(4))')*trace((d-eye(4)));
            
            v = v*d;% + G*W;
            v(:,1) = [];
            TY = v;
        end
        lambda1 = annealing*lambda1; 
        lambda2 = annealing*lambda2;
        sigma2 = sigma2*annealing;
        
        [az, el] = view;
        clf;
        plot3(X(:,1), X(:,2), X(:,3),'.r','MarkerSize',10);
        
        hold on;
        plot3(TY(:,1), TY(:,2), TY(:,3),'ob','MarkerSize', 5, 'MarkerFaceColor', 'b');
        legend({'target','source'},'location','NEo')
        view(az, el);
        drawnow;
        
        iters = iters+1;
    end

end

