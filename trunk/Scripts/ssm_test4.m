add_bcpd_paths;
clear; clc;

%% Load the surfaces
load('..\data\prostate\surfaces.mat');

%% SSM construction options
ssmOpt = struct;
ssmOpt.affine=1; ssmOpt.nonrigid=1; ssmOpt.beta=1;
ssmOpt.E=48.0; ssmOpt.nu=0.49; ssmOpt.save=1;
ssmOpt.viz=1;
ssmOpt.pointNr=300; ssmOpt.w=0.0;

% Stop conditions
ssmOpt.maxGiters=2; %Groupwise maximum number of iterations
ssmOpt.maxRiters=2; %Registration (rigid and non-rigid) maximum number of iterations
ssmOpt.sigmaTol=1e-6; %Sigma tolerance throughout registration
ssmOpt.QtTol=1e-6; %Error tolerance to update shape
ssmOpt.QgTol=1e-6; %Error tolerance to terminate groupwise registration

% Step sizes

[SSM correspondences] = cpd_ssmCon(fv, ssmOpt)