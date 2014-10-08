function R = rotateM_Jor(Rot)
%ROTATEM_JOR compute the rotation matrix given three Euler angles
%
% R = ROTATEM_JOR(Rot)
%
% Rot       (3 by 1)    Rotation angles
% 
% R         (3 by 3)    Rotation matrix
%
%   Written by Abtin Rasoulian (abtinr@ece.ubc.ca)
Rot = (Rot *  pi /180);
R =[cos(Rot(3)) -sin(Rot(3)) 0; sin(Rot(3)) cos(Rot(3)) 0; 0 0 1]*[cos(Rot(2)) 0 sin(Rot(2)); 0 1 0; -sin(Rot(2)) 0 cos(Rot(2))]*[1 0 0; 0 cos(Rot(1)) -sin(Rot(1)); 0 sin(Rot(1)) cos(Rot(1))];  




