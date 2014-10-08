add_bcpd_paths;

%% Load SSM
load 'SSM.mat';

% # modes to keep
K = 10;
filename = '../Data/prostate.ssm';

write_ssm(SSM.mean, SSM.mods(:,1:K), SSM.latent(1:K), SSM.faces, ...
    filename, '% .17g');
[SSM2.mean, SSM2.mods, SSM2.latent, SSM2.faces] = read_ssm(filename);

%% compare
meanErr = max(abs(SSM.mean(:)-SSM2.mean(:)));
disp('Mean error:');
disp(meanErr);

modsErr = max(reshape(abs(SSM.mods(:,1:K)-SSM2.mods), [], 1));
disp('Mods error:');
disp(modsErr);

latentErr = max(abs(SSM.latent(1:K)-SSM2.latent(:)));
disp('Latent error:');
disp(latentErr);

facesErr = max(reshape(abs(SSM.faces-SSM2.faces), [], 1));
disp('Faces error:');
disp(facesErr);