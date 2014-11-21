add_bcpd_paths;

%% load input data
clear; clc; close all;
load('..\data\SSM.mat');
[X,Fx] = readPolyDataInVTK('..\data\prostate\P1.vtk');

% [sphere.Vertices, sphere.Faces] = read_obj('sphere1.obj');
% [sphere2.Vertices, sphere2.Faces] = read_obj('sphere2.obj');

 clf;

% FV = sphere2;
% Y = sphere.Vertices;
% Fy = sphere.Faces;
 % X = sphere2.Vertices;
 % Fx = sphere2.Faces;
 
FV.faces = Fx;
FV.vertices = X;
Y = SSM.mean;
Fy = SSM.faces;



[az, el] = view;
clf;
patch(FV,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
hold on;
axis([-3,3,-3,3,-3,3]);
% without deformation
plot3(Y(:,1), Y(:,2), Y(:,3),'ob','MarkerSize', 5, 'MarkerFaceColor', 'b');
axis([-3,3,-3,3,-3,3]);
legend({'target','ssm'},'location','NEo')
view(az, el);
drawnow;

%% CPD values
errtol=1e-10; maxiters=50;
beta=1.0; E=48; nu = 1/3;
nMods=50; mu=3.0;

%% Do a rigid registration first
w = 0;
tic;
[TY, R, t, s, ~, sigma2 ] = cpd_rigid(X, Y, w, errtol, maxiters, eye(3), [0;0;0], 10.0, []);
% transform
Y = bsxfun(@plus, Y*(R')*s, t');
time = toc; fprintf('Elapsed rigid time is %f seconds\n', time);
fprintf('Rigid sigma is %f seconds\n', sigma2);

clf;
patch(FV,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
hold on;
axis([-3,3,-3,3,-3,3]);
% without deformation
plot3(Y(:,1), Y(:,2), Y(:,3),'ob','MarkerSize', 5, 'MarkerFaceColor', 'b');
axis([-3,3,-3,3,-3,3]);
legend({'target','ssm'},'location','NEo')
view(az, el);
drawnow;
%%
% Perform icp-fem registration between the two
beta=1e-1; E=1.0; nu=0.40; maxiters= 10000;

tic;
[TYd, ~, fem, u, ~] = icp_fem(X, TY, Fy, errtol, maxiters, eye(3), [0;0;0], 1.0, beta, E, nu, [], [], []);
time=toc; fprintf('FEM registration time is %f seconds\n', time);

%% Write out the registered SSM

% %% Calculate the dice after registration
% This section does not work, the meshes are not closed.
% V_fixed = abs(VolumeSurface(X,f));
% V_reg = abs(VolumeSurface(TY,SSM.faces));
% write_off('..\data\fixed.off',X,f); 
% write_off('..\data\reg.off',TY,SSM.faces);
% system('..\Core\VolumeOp.exe ..\data\fixed.off ..\data\reg.off ..\data\int_reg.off');
% [pint,fint] = read_off('..\data\int_reg.off');
% V_int_reg = abs(VolumeSurface(pint',fint'));
% Dice_after = 2.0*V_int_reg/(V_fixed+V_reg);
% 
% %% Calculate the dice before registration
% V_mov = abs(VolumeSurface(SSM.mean,SSM.faces));
% write_off('..\data\mov.off',SSM.mean,SSM.faces);
% system('..\Core\VolumeOp.exe ..\data\fixed.off ..\data\mov.off ..\data\int_mov.off');
% [pint,fint] = read_off('..\data\int_mov.off');
% V_int_mov = abs(VolumeSurface(pint',fint'));
% Dice_after = 2.0*V_int_mov/(V_fixed+V_mov);