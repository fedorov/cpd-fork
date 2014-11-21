add_bcpd_paths;

%% Read in the SSM and write out the mean shape as a VTK file
load('..\data\SSM.mat');
writePolyDataInVTK(SSM.mean, SSM.faces, '..\data\mean.vtk');

%% Read in the mean from the vtk file and compare it with SSM.mean
[v,f] = read_vtk('..\data\mean.vtk');

scatter3(SSM.mean(:,1), SSM.mean(:,2), SSM.mean(:,3), 'xr');
hold on;
scatter3(v(1,:)', v(2,:)', v(3,:)', 'ob');