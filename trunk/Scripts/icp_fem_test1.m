add_bcpd_paths;

clear; clc;

% Read the MR and US
X = read_fcsv('..\data\gmm-fem-trial2\P068_US_Seg.fcsv');
[Y,fY] = read_vrml('..\data\gmm-fem-trial2\P068_MR_Seg.wrl');

Y = Y'; fY = fY';

R = rotation_matrix( 0, 90, 180, 'degrees' );

Y = Y*R';

% Do a center of mass alignment
t = (mean(X) - mean(Y))';
Y = bsxfun(@plus,Y,t');

% Perform a rigid registration between the two
w=0.00; errtol=1e-4; maxiters=100;

[TY, R, t, s, ~, ~ ] = cpd_rigid(X, Y, w, errtol, maxiters, eye(3), [0;0;0], 10.0, []);

% Perform icp-fem registration between the two
beta=1e-1; E=1.0; nu=0.40; maxiters= 10000;

tic;
[TYd, ~, fem, u, ~] = icp_fem(X, TY, fY, errtol, maxiters, eye(3), [0;0;0], 1.0, beta, E, nu, [], [], []);
time=toc; fprintf('FEM registration time is %f seconds\n', time);