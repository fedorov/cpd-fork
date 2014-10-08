function [b, sigma2] = cpd_ssm( SSM, X, f, numModes, b, w, errtol, maxiters, lambda, sigma2)
%CPD_SSM finds the correct SSM modes that warp two point sets together
%   b = cpd_ssm( SSM, X, b, numModes, w, errtol, maxiters, sigma2)
%
%   Inputs:
%   SSM  Structure to hold the SSM
%   X    M x D matrix of data points
%   f    M x D matrix of faces
%   numModes    Number of Modes to use in registration
%   b    Weights that control the shape of the SSM
%   w    Weight to account for noise and outliers (optional, defaults to 0)
%   errtol Error tolerance for convergence.  Algorithm will terminate if
%        the objective function does not change by more than this amount
%        (optional, defaults to 1e-10)
%   maxiters Maximum number of iterations.  Algorithm will terminate if
%        this many iterations have been performed (optional, defaults to
%        100)
%   sigma2  Initial estimate of squared variance (optional, default
%        estimated from SSM, X)
%   
%
%   Outputs:

[M,D] = size(SSM.mean);
N = size(X,1);

if (nargin < 5 || isempty(b))
    b = zeros(1, numModes);
end
if (nargin < 6 || isempty(w))
    w=0.0;
end
if (nargin < 7 || isempty(errtol))
    errtol = 1e-10;
end

if (nargin < 8 || isempty(maxiters))
    maxiters=100;
end

%Initialize the instance
ch = SSM.mods(:,1:numModes)*b';
ad = reshape(ch, 3, length(ch)/3)';
Y = SSM.mean + ad;

if(nargin < 8 || isempty(lambda))
    lambda = 0.0;
end

if (nargin < 9 || isempty(sigma2))
    XX = reshape(X, [1, N, D]);
    YY = reshape(Y, [M, 1, D]);
    XX = repmat(XX, [M, 1, 1]);
    YY = repmat(YY, [1, N, 1]);
    diff = XX-YY;
    diff = diff.*diff;
    err = sum(diff(:));
    sigma2 = 1/(D*N*M)*err;
    clear diff err;
end

m_mods = SSM.mods(:,1:numModes);
m_latent = SSM.latent(1:numModes);

iter = 0;

while( iter < maxiters && sigma2 > errtol )
    %E-step
    [P,P1,Pt1,~] = cpd_P( X, Y, sigma2, w);
    
    %M-step
    b = findOptimumModsMult(Y, X, P, m_mods, m_latent,lambda);
    ch = m_mods*b;
    ad = reshape(ch, 3, length(ch)/3)';
    Y = Y+ad;
    
    %Update sigma
    XX = reshape(X, [1, N, D]);
    YY = reshape(Y, [M, 1, D]);
    XX = repmat(XX, [M, 1, 1]);
    YY = repmat(YY, [1, N, 1]);
    diff = XX-YY;
    diff = diff.*diff;
    err = sum(diff(:));
    sigma2 = 1/(D*N*M)*err;
    clear diff err;
    iter = iter + 1;
end

end

