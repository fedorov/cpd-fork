% Generates different shape instances of the prostate, and regenerates a
% tet mesh if required.

add_bcpd_paths

%% Load prostate SSM and generate a model

% modes to keep
K = 10;

load SSM.mat;
X = reshape(SSM.mean',[],1);
P = SSM.mods(:,1:K);
S = SSM.latent(1:K);
F = SSM.faces;

% Tet model of mean
Y = SSM.mean;
[N, ~, E] = tetgen_mex(Y', F',[],'');
prostate = fem_model(N', E');

D = fem_material_linear(44000, 0.4);

% Tell prostate to abort stiffness computation for inverted elements
setAbortJacobian(prostate, 1e-10);    

%% Loop generating prostate fem

Niters = 50;
for i=1:Niters
    % generate new shape instance
    b = 0.5*rand(K,1).*S;
    Y = reshape(X+P*b,3,[])';
    
    % update nodes in existing model
    updateNodes(prostate,Y, 1:size(Y,1));
    % compute stiffness
    [Km, minJ] = getStiffnessMatrix(prostate, D);
    if (isempty(Km))
        fprintf('Inverted element detected, regenerating mesh... ');
        tic;
        [N, ~, E] = tetgen_mex(Y', F', 2, '');
        rtime = toc;
        fprintf('%g s\n', rtime);
        set(prostate, N', E');
        [Km, minJ] = getStiffnessMatrix(prostate, D);
    else
        fprintf('Mesh can be re-used, minJ=%f\n',minJ);
    end
    
    hold off;
    N = getNodes(prostate);
    plot3(N(:,1), N(:,2), N(:,3),'.k','MarkerSize',10);
    hold on;
    pmesh.vertices = Y;
    pmesh.faces = F;
    patch(pmesh,'FaceColor','red','FaceAlpha',0.1);
    axis([-2 2 -2 2 -2 2]);
    hold off;
    
    pause(0.1);
    
end
