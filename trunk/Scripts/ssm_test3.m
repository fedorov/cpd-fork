add_bcpd_paths;
clear; clc;

%% Read in the SSM
load('..\data\SSM.mat');

%% Read in the segmented prostate
[X,f] = readPolyDataInVTK('..\data\prostate\P1.vtk');
Y = SSM.mean;

A=eye(3); t=[0;0;0];

%% Affine registration
for i=1:10
    i
    if (i==1)
        sigma2=[];
    end
    tic;
    [~, A, t, ~, sigma2] = cpd_affine(X, Y, [], 1e-5, 100, A,t,[]);
    time = toc; fprintf('Elapsed affine time is %f seconds\n', time);
    fprintf('Affine sigma is %f seconds\n', sigma2);
    
    %% Write out the affine registration results
%     Y = bsxfun(@plus, SSM.mean*(Aout'), tout');
%     writePolyDataInVTK(Y, SSM.faces,'..\data\P1_aff.vtk');
    
    %% Correct the fixed image with the affine transform and write it out
    Xb = bsxfun(@minus, X*(inv(A')), t');
%     writePolyDataInVTK(Xb, f,'..\data\P1_b.vtk');
    
    %% Perform SSM registration
    numModes = 50;
    tic;
    [b,sigma2] = cpd_ssm( SSM, Xb, f, numModes, [], [], [], 100, 0, sigma2);
    time = toc; fprintf('Elapsed shape time is %f seconds\n', time);
    fprintf('SSM sigma is %f seconds\n', sigma2);
    
    %% Write out the shape registration results
    ch = SSM.mods(:,1:numModes)*b;
    ad = reshape(ch, 3, length(ch)/3)';
    Y = SSM.mean+ad;
    Ytemp = bsxfun(@plus, Y*(A'), t');
    writePolyDataInVTK(Ytemp, SSM.faces,'..\data\P1_temp_ssm.vtk');
end
writePolyDataInVTK(Y, SSM.faces,'..\data\P1_ssm.vtk');