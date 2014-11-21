add_bcpd_paths;

%% Read in the SSM
load('..\data\SSM.mat');

%% Test properties
N = 5; %number of modes to use
w = [1.0 0.5 0.25 0.3 0.5]';
w = SSM.latent(1:N).*w; %latent holds the eigenvalues of the SSM in descending order

%% Create a new instance
ch = SSM.mods(:,1:N)*w;
ad = reshape(ch, 3, length(ch)/3)';
I = SSM.mean + ad;

%% Plot the new instance of the SSM against the mean
clf;
scatter3(I(:,1), I(:,2), I(:,3), 'xr');
hold on;
scatter3(SSM.mean(:,1), SSM.mean(:,2), SSM.mean(:,3), 'ob');

%% Call atlas viewer to do create and visualize instances interactively
atlasViewer(SSM); % Press the "Create Plot" button to visualize