function R = cpd_Rotation( alpha, beta, gamma, mode )
%CPD_ROTATION Creates a 3 x 3 rotation matrix
%   R = cpd_Rotation( alpha, beta, gamma, mode )
%
%   Inputs:
%   alpha  Rotation around X-axis, default radians
%   beta   Rotation around Y-axis, default radians
%   gamma  Rotation around Z-axis, default radians
%   mode   Rotations expressed in degrees or radians

%   Outputs:
%   R      3 x 3 rotation matrix

if (nargin < 4 || strcmp(mode, '') || isempty(mode))
    mode = 'radians';
end

if ( isempty(alpha) );
    alpha = 0.0;
end

if ( isempty(beta) );
    beta = 0.0;
end

if ( isempty(gamma) );
    gamma = 0.0;
end

if ( strcmp(mode, 'degrees') )
    alpha = alpha*pi/180.0; beta = beta*pi/180.0; gamma = gamma*pi/180.0;
end

Rx = eye(3); Ry = eye(3); Rz = eye(3);

Rx(2,2) = cos(alpha); Rx(2,3) = -sin(alpha); Rx(3,2) = sin(alpha); Rx(3,3)= cos(alpha);
Ry(1,1) = cos(beta); Ry(1,3) = sin(beta); Ry(3,1) = -sin(beta); Ry(3,3)= cos(beta);
Rz(1,1) = cos(gamma); Rz(1,2) = -sin(gamma); Rz(2,1) = sin(gamma); Rz(2,2)= cos(gamma);

R = Rx*Ry*Rz;

end

