add_bcpd_paths;

clear; clc; close all;

% store the errors
TRE_before = [];
TRE_after = [];

Dice_before = [];
Dice_after = [];

%% Read in the label maps
[Y,fY] = read_ply('C:\data\ProstateQueens\P030-data\US_affine_model.ply');
[X,fX] = read_ply('C:\data\ProstateQueens\P030-data\MR_model.ply');

numberOfFaces = 2000;

pX.faces = fX;
pX.vertices = X;

pY.faces = fY;
pY.vertices = Y;

landmark_MR = [4.64 -70.300 -1.825; 4.579 -68.214 -3.795; -2.736 -58.524 -7.623; -1.984 -60.930 -7.623];
landmark_US = [5.281 -69.204 -1.825; 4.697 -66.087 -3.795; -3.112 -56.945 -7.623; -2.510 -58.975 -7.623];

% [fX,X] = reducepatch(pX,numberOfFaces);
nfx = reducepatch(pX,numberOfFaces);
nfy = reducepatch(pY,numberOfFaces);
figure, patch(nfx,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
hold on;
patch(nfy,'FaceColor','blue','FaceAlpha',0.2, 'EdgeColor', 'blue', 'EdgeAlpha',0.2);
hold on;
plot3(landmark_MR(:,1), landmark_MR(:,2), landmark_MR(:,3),'og','MarkerSize', 5, 'MarkerFaceColor', 'g');
plot3(landmark_US(:,1), landmark_US(:,2), landmark_US(:,3),'xg','MarkerSize', 5, 'MarkerFaceColor', 'g');

%% Calculate dice using the "inpolyhedron" technique
tic;
dice_polyhedron = cpd_dice(nfx.vertices,nfx.faces,nfy.vertices,nfy.faces,50,50,50);
dice_time = toc;

Dice_before = [Dice_before; dice_polyhedron];

%% Do the registration
w=0.00; errtol=1e-4; maxiters=500; sigma2=10;
beta=0.12; E=4.8; nu=0.49;

tic;
[TY, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(nfx.vertices, nfy.vertices, nfy.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu, [], [], []);
time=toc; fprintf('FEM registration time is %f seconds\n', time);

tic;
Phi = getInterpolationMatrix(fem, landmark_US);
time=toc; fprintf('FEM interpolation time is %f seconds\n', time);

landmark_warped_US = landmark_US + Phi*u;

error_before_reg = sqrt(sum((landmark_MR - landmark_US).*(landmark_MR - landmark_US),2));
error_after_reg = sqrt(sum((landmark_MR - landmark_warped_US).*(landmark_MR - landmark_warped_US),2));

TRE_before = [TRE_before; error_before_reg];
TRE_after = [TRE_after; error_after_reg];

tic;
dice_polyhedron = cpd_dice(nfx.vertices,nfx.faces,TY,nfy.faces,50,50,50);
dice_time = toc;

Dice_after = [Dice_after; dice_polyhedron];
%% Read in the label maps
[Y,fY] = read_ply('C:\data\ProstateQueens\P032-data\US_affine_model.ply');
[X,fX] = read_ply('C:\data\ProstateQueens\P032-data\MR_model.ply');

numberOfFaces = 2000;

pX.faces = fX;
pX.vertices = X;

pY.faces = fY;
pY.vertices = Y;

landmark_MR = [-1.025 -56.746 38.967];
landmark_US = [-0.932 -58.874 38.967];

nfx = reducepatch(pX,numberOfFaces);
nfy = reducepatch(pY,numberOfFaces);
figure, patch(nfx,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
hold on;
patch(nfy,'FaceColor','blue','FaceAlpha',0.2, 'EdgeColor', 'blue', 'EdgeAlpha',0.2);
hold on;
plot3(landmark_MR(:,1), landmark_MR(:,2), landmark_MR(:,3),'og','MarkerSize', 5, 'MarkerFaceColor', 'g');
plot3(landmark_US(:,1), landmark_US(:,2), landmark_US(:,3),'xg','MarkerSize', 5, 'MarkerFaceColor', 'g');

tic;
dice_polyhedron = cpd_dice(nfx.vertices,nfx.faces,nfy.vertices,nfy.faces,50,50,50);
dice_time = toc;

Dice_before = [Dice_before; dice_polyhedron];

%% Do the registration
w=0.00; errtol=1e-4; maxiters=500; sigma2=10;
beta=0.01; E=4.8; nu=0.49;

tic;
[TY, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(nfx.vertices, nfy.vertices, nfy.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu, [], [], []);
time=toc; fprintf('FEM registration time is %f seconds\n', time);

tic;
Phi = getInterpolationMatrix(fem, landmark_US);
time=toc; fprintf('FEM interpolation time is %f seconds\n', time);

landmark_warped_US = landmark_US + Phi*u;

error_before_reg = sqrt(sum((landmark_MR - landmark_US).*(landmark_MR - landmark_US),2));
error_after_reg = sqrt(sum((landmark_MR - landmark_warped_US).*(landmark_MR - landmark_warped_US),2));

TRE_before = [TRE_before; error_before_reg];
TRE_after = [TRE_after; error_after_reg];

tic;
dice_polyhedron = cpd_dice(nfx.vertices,nfx.faces,TY,nfy.faces,50,50,50);
dice_time = toc;

Dice_after = [Dice_after; dice_polyhedron];

%% Read in the label maps
[Y,fY] = read_ply('C:\data\ProstateQueens\P034-data\US_affine_model.ply');
[X,fX] = read_ply('C:\data\ProstateQueens\P034-data\MR_model.ply');

numberOfFaces = 1000;

pX.faces = fX;
pX.vertices = X;

pY.faces = fY;
pY.vertices = Y;

landmark_MR = [6.666 -52.673 -2.195; 9.478 -56.695 7.237; 11.578 -57.084 7.237; 7.922 -56.228 5.323; 9.555 -56.850 5.323; 4.265 -56.073 3.136];
landmark_US = [8.683 -51.762 -2.195; 10.628 -52.680 7.237; 13.024 -54.096 7.237; 12.806 -53.334 5.323; 15.638 -53.334 5.323; 7.144 -54.423 3.674];

nfx = reducepatch(pX,numberOfFaces);
nfy = reducepatch(pY,numberOfFaces);
figure, patch(nfx,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
hold on;
patch(nfy,'FaceColor','blue','FaceAlpha',0.2, 'EdgeColor', 'blue', 'EdgeAlpha',0.2);
hold on;
plot3(landmark_MR(:,1), landmark_MR(:,2), landmark_MR(:,3),'og','MarkerSize', 5, 'MarkerFaceColor', 'g');
plot3(landmark_US(:,1), landmark_US(:,2), landmark_US(:,3),'xg','MarkerSize', 5, 'MarkerFaceColor', 'g');

tic;
dice_polyhedron = cpd_dice(nfx.vertices,nfx.faces,nfy.vertices,nfy.faces,50,50,50);
dice_time = toc;

Dice_before = [Dice_before; dice_polyhedron];

%% Do the registration
w=0.00; errtol=1e-4; maxiters=500; sigma2=10;
beta=0.01; E=4.8; nu=0.49;

tic;
[TY, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(nfx.vertices, nfy.vertices, nfy.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu, [], [], []);
time=toc; fprintf('FEM registration time is %f seconds\n', time);

tic;
Phi = getInterpolationMatrix(fem, landmark_US);
time=toc; fprintf('FEM interpolation time is %f seconds\n', time);

landmark_warped_US = landmark_US + Phi*u;

error_before_reg = sqrt(sum((landmark_MR - landmark_US).*(landmark_MR - landmark_US),2));
error_after_reg = sqrt(sum((landmark_MR - landmark_warped_US).*(landmark_MR - landmark_warped_US),2));

TRE_before = [TRE_before; error_before_reg];
TRE_after = [TRE_after; error_after_reg];

tic;
dice_polyhedron = cpd_dice(nfx.vertices,nfx.faces,TY,nfy.faces,50,50,50);
dice_time = toc;

Dice_after = [Dice_after; dice_polyhedron];

%% Read in the label maps
[Y,fY] = read_ply('C:\data\ProstateQueens\P035-data\US_affine_model.ply');
[X,fX] = read_ply('C:\data\ProstateQueens\P035-data\MR_model.ply');

numberOfFaces = 2000;

pX.faces = fX;
pX.vertices = X;

pY.faces = fY;
pY.vertices = Y;

landmark_MR = [7.714 -92.983 3.326; 3.972 -91.927 3.326; 6.275 -88.377 3.326; 8.481 -81.852 4.693];
landmark_US = [3.684 -91.351 3.326; 0.422 -89.049 3.326; 5.027 -86.458 3.326; 7.618 -80.605 4.346];

nfx = reducepatch(pX,numberOfFaces);
nfy = reducepatch(pY,numberOfFaces);
figure, patch(nfx,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
hold on;
patch(nfy,'FaceColor','blue','FaceAlpha',0.2, 'EdgeColor', 'blue', 'EdgeAlpha',0.2);
hold on;
plot3(landmark_MR(:,1), landmark_MR(:,2), landmark_MR(:,3),'og','MarkerSize', 5, 'MarkerFaceColor', 'g');
plot3(landmark_US(:,1), landmark_US(:,2), landmark_US(:,3),'xg','MarkerSize', 5, 'MarkerFaceColor', 'g');

tic;
dice_polyhedron = cpd_dice(nfx.vertices,nfx.faces,nfy.vertices,nfy.faces,50,50,50);
dice_time = toc;

Dice_before = [Dice_before; dice_polyhedron];

%% Do the registration
w=0.00; errtol=1e-4; maxiters=500; sigma2=10;
beta=0.01; E=4.8; nu=0.49;

tic;
[TY, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(nfx.vertices, nfy.vertices, nfy.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu, [], [], []);
time=toc; fprintf('FEM registration time is %f seconds\n', time);

tic;
Phi = getInterpolationMatrix(fem, landmark_US);
time=toc; fprintf('FEM interpolation time is %f seconds\n', time);

landmark_warped_US = landmark_US + Phi*u;

error_before_reg = sqrt(sum((landmark_MR - landmark_US).*(landmark_MR - landmark_US),2));
error_after_reg = sqrt(sum((landmark_MR - landmark_warped_US).*(landmark_MR - landmark_warped_US),2));

TRE_before = [TRE_before; error_before_reg];
TRE_after = [TRE_after; error_after_reg];

tic;
dice_polyhedron = cpd_dice(nfx.vertices,nfx.faces,TY,nfy.faces,50,50,50);
dice_time = toc;

Dice_after = [Dice_after; dice_polyhedron];

% %% Read in the label maps
% [Y,fY] = read_ply('C:\data\ProstateQueens\P036-data\US_affine_model.ply');
% [X,fX] = read_ply('C:\data\ProstateQueens\P036-data\MR_model.ply');
% 
% numberOfFaces = 2000;
% 
% pX.faces = fX;
% pX.vertices = X;
% 
% pY.faces = fY;
% pY.vertices = Y;
% 
% landmark_MR = [4.150 -85.329 4.633];
% landmark_US = [3.891 -85.329 3.725];
% 
% nfx = reducepatch(pX,numberOfFaces);
% nfy = reducepatch(pY,numberOfFaces);
% figure, patch(nfx,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
% hold on;
% patch(nfy,'FaceColor','blue','FaceAlpha',0.2, 'EdgeColor', 'blue', 'EdgeAlpha',0.2);
% hold on;
% 
% %% Do the registration
% w=0.00; errtol=1e-4; maxiters=500; sigma2=10;
% beta=0.01; E=4.8; nu=0.49;
% 
% tic;
% [TY, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(nfx.vertices, nfy.vertices, nfy.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu, [], [], []);
% time=toc; fprintf('FEM registration time is %f seconds\n', time);
% 
% tic;
% Phi = getInterpolationMatrix(fem, landmark_US);
% time=toc; fprintf('FEM interpolation time is %f seconds\n', time);
% 
% landmark_warped_US = landmark_US + Phi*u;
% 
% error_before_reg = sqrt(sum((landmark_MR - landmark_US).*(landmark_MR - landmark_US),2));
% error_after_reg = sqrt(sum((landmark_MR - landmark_warped_US).*(landmark_MR - landmark_warped_US),2));
% 
% TRE_before = [TRE_before; error_before_reg];
% TRE_after = [TRE_after; error_after_reg];

%% Read in the label maps
[Y,fY] = read_ply('C:\data\ProstateQueens\P2002-data\US_affine_model.ply');
[X,fX] = read_ply('C:\data\ProstateQueens\P2002-data\MR_model.ply');

numberOfFaces = 2000;

pX.faces = fX;
pX.vertices = X;

pY.faces = fY;
pY.vertices = Y;

landmark_MR = [7.942 -80.861 11.278; 6.442 -79.534 -12.470];
landmark_US = [7.042 -78.224 13.707; 8.377 -77.406 -11.544];

nfx = reducepatch(pX,numberOfFaces);
nfy = reducepatch(pY,numberOfFaces);
figure, patch(nfx,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
hold on;
patch(nfy,'FaceColor','blue','FaceAlpha',0.2, 'EdgeColor', 'blue', 'EdgeAlpha',0.2);
hold on;
plot3(landmark_MR(:,1), landmark_MR(:,2), landmark_MR(:,3),'og','MarkerSize', 5, 'MarkerFaceColor', 'g');
plot3(landmark_US(:,1), landmark_US(:,2), landmark_US(:,3),'xg','MarkerSize', 5, 'MarkerFaceColor', 'g');


tic;
dice_polyhedron = cpd_dice(nfx.vertices,nfx.faces,nfy.vertices,nfy.faces,50,50,50);
dice_time = toc;

Dice_before = [Dice_before; dice_polyhedron];

%% Do the registration
w=0.00; errtol=1e-4; maxiters=500; sigma2=10;
beta=0.01; E=4.8; nu=0.49;

tic;
[TY, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(nfx.vertices, nfy.vertices, nfy.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu, [], [], []);
time=toc; fprintf('FEM registration time is %f seconds\n', time);

tic;
Phi = getInterpolationMatrix(fem, landmark_US);
time=toc; fprintf('FEM interpolation time is %f seconds\n', time);

landmark_warped_US = landmark_US + Phi*u;

error_before_reg = sqrt(sum((landmark_MR - landmark_US).*(landmark_MR - landmark_US),2));
error_after_reg = sqrt(sum((landmark_MR - landmark_warped_US).*(landmark_MR - landmark_warped_US),2));

TRE_before = [TRE_before; error_before_reg];
TRE_after = [TRE_after; error_after_reg];

tic;
dice_polyhedron = cpd_dice(nfx.vertices,nfx.faces,TY,nfy.faces,50,50,50);
dice_time = toc;

Dice_after = [Dice_after; dice_polyhedron];

%% Read in the label maps
[Y,fY] = read_ply('C:\data\ProstateQueens\P2004-data\US_affine_model.ply');
[X,fX] = read_ply('C:\data\ProstateQueens\P2004-data\MR_model.ply');

numberOfFaces = 1500;

pX.faces = fX;
pX.vertices = X;

pY.faces = fY;
pY.vertices = Y;

landmark_MR = [2.162 -80.535 -5.629; 2.186 -80.116 -5.778; 2.186 -78.298 -12.271];
landmark_US = [2.162 -76.426 -5.629; 2.186 -78.125 -8.548; 2.186 -74.316 -14.868];

nfx = reducepatch(pX,numberOfFaces);
nfy = reducepatch(pY,numberOfFaces);
figure, patch(nfx,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
hold on;
patch(nfy,'FaceColor','blue','FaceAlpha',0.2, 'EdgeColor', 'blue', 'EdgeAlpha',0.2);
hold on;
% plot3(landmark_MR(:,1), landmark_MR(:,2), landmark_MR(:,3),'og','MarkerSize', 5, 'MarkerFaceColor', 'g');
% plot3(landmark_US(:,1), landmark_US(:,2), landmark_US(:,3),'xg','MarkerSize', 5, 'MarkerFaceColor', 'g');

tic;
dice_polyhedron = cpd_dice(nfx.vertices,nfx.faces,nfy.vertices,nfy.faces,50,50,50);
dice_time = toc;

Dice_before = [Dice_before; dice_polyhedron];

%% Do the registration
w=0.00; errtol=1e-4; maxiters=500; sigma2=10;
beta=0.01; E=4.8; nu=0.49;

tic;
[TY, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(nfx.vertices, nfy.vertices, nfy.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu, [], [], []);
time=toc; fprintf('FEM registration time is %f seconds\n', time);

tic;
Phi = getInterpolationMatrix(fem, landmark_US);
time=toc; fprintf('FEM interpolation time is %f seconds\n', time);

landmark_warped_US = landmark_US + Phi*u;

error_before_reg = sqrt(sum((landmark_MR - landmark_US).*(landmark_MR - landmark_US),2));
error_after_reg = sqrt(sum((landmark_MR - landmark_warped_US).*(landmark_MR - landmark_warped_US),2));

TRE_before = [TRE_before; error_before_reg];
TRE_after = [TRE_after; error_after_reg];

tic;
dice_polyhedron = cpd_dice(nfx.vertices,nfx.faces,TY,nfy.faces,50,50,50);
dice_time = toc;

Dice_after = [Dice_after; dice_polyhedron];


%% Read in the label maps
[Y,fY] = read_ply('C:\data\ProstateQueens\P029-data\US_affine_model.ply');
[X,fX] = read_ply('C:\data\ProstateQueens\P029-data\MR_model.ply');

numberOfFaces = 2000;

pX.faces = fX;
pX.vertices = X;

pY.faces = fY;
pY.vertices = Y;

landmark_MR = [-4.455 -49.867 -4.168];
landmark_US = [-3.449 -51.647 -4.168];

nfx = reducepatch(pX,numberOfFaces);
nfy = reducepatch(pY,numberOfFaces);
figure, patch(nfx,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
hold on;
patch(nfy,'FaceColor','blue','FaceAlpha',0.2, 'EdgeColor', 'blue', 'EdgeAlpha',0.2);
hold on;
% plot3(landmark_MR(:,1), landmark_MR(:,2), landmark_MR(:,3),'og','MarkerSize', 5, 'MarkerFaceColor', 'g');
% plot3(landmark_US(:,1), landmark_US(:,2), landmark_US(:,3),'xg','MarkerSize', 5, 'MarkerFaceColor', 'g');

tic;
dice_polyhedron = cpd_dice(nfx.vertices,nfx.faces,nfy.vertices,nfy.faces,50,50,50);
dice_time = toc;

Dice_before = [Dice_before; dice_polyhedron];

%% Do the registration
w=0.00; errtol=1e-4; maxiters=500; sigma2=10;
beta=0.01; E=4.8; nu=0.49;

tic;
[TY, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(nfx.vertices, nfy.vertices, nfy.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu, [], [], []);
time=toc; fprintf('FEM registration time is %f seconds\n', time);

tic;
Phi = getInterpolationMatrix(fem, landmark_US);
time=toc; fprintf('FEM interpolation time is %f seconds\n', time);

landmark_warped_US = landmark_US + Phi*u;

error_before_reg = sqrt(sum((landmark_MR - landmark_US).*(landmark_MR - landmark_US),2));
error_after_reg = sqrt(sum((landmark_MR - landmark_warped_US).*(landmark_MR - landmark_warped_US),2));

TRE_before = [TRE_before; error_before_reg];
TRE_after = [TRE_after; error_after_reg];

tic;
dice_polyhedron = cpd_dice(nfx.vertices,nfx.faces,TY,nfy.faces,50,50,50);
dice_time = toc;

Dice_after = [Dice_after; dice_polyhedron];