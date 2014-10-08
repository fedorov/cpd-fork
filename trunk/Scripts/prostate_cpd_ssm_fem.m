add_bcpd_paths;

%% load input data
clear; clc; close all;
load('..\data\SSM.mat');
[X,f] = readPolyDataInVTK('..\data\prostate\P1.vtk');

 clf;
FV.faces = f;
FV.vertices = X;

[az, el] = view;
clf;
patch(FV,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
hold on;
axis([-3,3,-3,3,-3,3]);
% without deformation
plot3(X(:,1), X(:,2), X(:,3),'ob','MarkerSize', 5, 'MarkerFaceColor', 'b');
axis([-3,3,-3,3,-3,3]);
legend({'target','ssm'},'location','NEo')
view(az, el);
drawnow;

%% CPD values
s = 1; R=eye(3); t=[0;0;0];
w=0.0; errtol=1e-10; maxiters=50; sigma2=[];
beta=1.0; E=480; nu=1/3;
nMods=50; mu=3.0;

%% Do a rigid registration first
% tic;
% [~, R, t, s, ~, sigma2] = cpd_rigid(X, SSM.mean, w, errtol, maxiters, [],[],[], sigma2);
% time = toc; fprintf('Elapsed rigid time is %f seconds\n', time);
% fprintf('Rigid sigma is %f seconds\n', sigma2);

%% Perform Biomechanically Constrained Point Drift (BCPD) registration
tic;
[TY, R, t, s, V, sigma2, P] = cpd_ssm_fem(X, SSM, w, errtol, maxiters, R, t, s, nMods, [], sigma2, mu, beta, E, nu, FV );
time=toc; fprintf('Registration time is %f seconds\n', time);
fprintf('Registration sigma2 is %f\n', sigma2);

%% Write out the registered SSM
TY = bsxfun(@plus, TY*(R')*s, t');
writePolyDataInVTK(TY,SSM.faces,'..\data\SSM_after.vtk');
writePolyDataInVTK(SSM.mean,SSM.faces,'..\data\SSM_before.vtk');

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