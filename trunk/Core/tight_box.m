function [ R, c, w, B ] = tight_box( X )
%TIGHT_BOX Computes a tight-fitting box for a collection of points
%   [R, c, w, B] = tight_box( X )
%       Computes a tight-fitting box around the collection of points X.
%       Uses an SVD to determine dominant axes, projects points along these
%       axes to determine limits.
%       
%       X:  Nx3 set of points
%       R:  Rotation matrix giving transformed box axes, defined by
%           post-multiplication
%       c:  center of the box
%       w:  widths of the box
%       B:  8x3 corners of the box, arranged with bottom face first, cw
%           w.r.t. outward normal, followed by corresponding corners
%           on top face
%
%       A box, B, with corners [+/- w(1)/2, +/- w(2)/2, +/- w(3)/2] will 
%       bound all points when transformed:
%           X in B*R+c

    c = mean(X,1);
    Y = bsxfun(@minus, X, c);
    D = (Y')*Y;
    [R, ~, ~] = svd(D);
    
    % Flip axis to be consistent with rotation
    if (det(R) < 0)
        R(:,3) = -R(:,3);
    end
    R = R'; % for post-multiply
    
    % Rotate to find projections
    Y = Y*(R');
    maxY = max(Y,[],1);
    minY = min(Y,[],1);
    w = (maxY-minY);
    
    % update centre of box
    c = (maxY+minY)/2*(R) + c;
    
    % compute corners of box
    B = [-1 -1 -1;
          1 -1 -1;
          1  1 -1;
         -1  1 -1;
         -1 -1  1;
          1 -1  1;
          1  1  1;
         -1  1  1 ];
     B = B.*repmat(w/2,[8 1]);
     B = B*R + repmat(c, [8, 1]);
               
end

