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
X = SSM.mean + ad;

%% Register the SSM to the instance
b = cpd_ssm( SSM, X, N, [], [], [], 100, [], []);