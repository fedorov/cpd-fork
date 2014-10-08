function [ R, t, s, TY ] = svd_reg( X, Y, mode )
%TIGHT_BOX_REG Registration of tight-fitting boxes
%   [R, t, s, TY] = tight_box_reg(X, Y, mode)
%       Registers the points Y to those in X by first computing a
%       tight-fitting box around each point set, and then registering those
%       boxes
%   X:    Nx3 reference points
%   Y:    Mx3 points to register
%   mode: 'rigid' or 'affine', determines whether the returned scale is a
%         scalar or a 3x1 vector (optional, default = 'rigid')
%   R:    Best-fitting rotation, post-multiplication
%   t:    translation
%   s:    scale, either scalar (rigid) or 3x1 vector (affine)
%   TY:   transformed points

    D = size(X,2);
    N = size(X,1);
    M = size(Y,1);

    modeIdx = 0;    % rigid
    if (nargin >= 3  && ~isempty(mode))
        if (mode(1) == 'a')
            modeIdx = 1;
        end
    end

    c1 = mean(X,1);
    XX = bsxfun(@minus, X, c1);
    [R1, w1, ~] = svd((XX')*XX);
    w1 = diag(w1)';
    % Flip axis to be consistent with rotation
    if (det(R1) < 0)
        R1(:,3) = -R1(:,3);
    end
    R1 = R1'; % for post-multiply
    
    c2 = mean(Y,1);
    XX = bsxfun(@minus, Y, c2);
    [R2, w2, ~] = svd((XX')*XX);
    w2 = diag(w2)';
    % Flip axis to be consistent with rotation
    if (det(R2) < 0)
        R2(:,3) = -R2(:,3);
    end
    R2 = R2'; % for post-multiply
    
    if (modeIdx == 0)
        % rigid scale
        s = sqrt(sum(w1)/sum(w2));
    else
        % affine scale
        s = sqrt(w1./w2)';
    end
   
    % default rotation
    R = (R2')*R1;
   
     % 4 possible rotations
    scales = [1 1 1; -1 1 -1; -1 -1 1; 1 -1 -1];
    minErr = Inf;
    for i=1:size(scales,1);
        % rotate axes
        Rtmp = diag(scales(i,:));
        Rtmp = (R2')*Rtmp*R1;    
        
        % transform
        if (modeIdx == 0)
            t = c1 - c2*Rtmp*s;
            TY = bsxfun(@plus, s*(Y)*(Rtmp), t);
        else
            t = c1 - c2*Rtmp*diag(s);
            TY = bsxfun(@plus, (Y)*(Rtmp)*diag(s), t);
        end
        
        XX = reshape(X, [1, N, D]);
        YY = reshape(TY, [M, 1, D]);
        XX = repmat(XX, [M, 1, 1]);
        YY = repmat(YY, [1, N, 1]);
        diff = XX-YY;
        diff = diff.*diff;
        err = sum(diff(:));
        
        if (err < minErr)
            R = Rtmp;
            minErr = err;
        end
        
    end

    % transform
    if (modeIdx == 0)
        t = c1 - c2*R*s;
        TY = bsxfun(@plus, s*(Y)*(R), t);
    else
        t = c1 - c2*R*diag(s);
        TY = bsxfun(@plus, (Y)*(R)*diag(s), t);
    end

end

