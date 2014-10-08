add_bcpd_paths;

%% Load and modify data
testcase = 'P30';

[X,Xf] = read_ply([testcase, '\MR_model.ply']);
[Y,Yf] = read_ply([testcase, '\US_affine_model.ply']);
maxFaces = 2000;

pX.faces = Xf;
pX.vertices = X;

pY.faces = Yf;
pY.vertices = Y;

nfx = reducepatch(pX,maxFaces);
nfy = reducepatch(pY,maxFaces);

X = nfx.vertices;
Xf = nfx.faces;
Y = nfy.vertices;
Yf = nfy.faces;

cd(testcase);
landmarks;
cd '..'

%% CPD-FEM
w=0.001; errtol=1e-4; maxiters=50; sigma2=10;
beta=0.12; E=4.8; nu=0.49; 

tic;
[TY, V, results(1).sigma2, ~, fem, U, Phi] = cpd_fem_norigid(X, Y, Yf, w, errtol, maxiters, sigma2, beta, E, nu);
% [TY, ~, ~, ~, ~, sigma2, ~, fem, U, ~] = cpd_fem_only(X, Y, Yf, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu, [], [], []); 
time=toc; fprintf('Registration time is %f seconds\n', time);
fprintf('Registration sigma2 is %f\n', results(1).sigma2);


%% Landmark error

% original FEM
setSurfaceMesh(fem, Yf);

% new FEM with moved node locations (for reverse interpolation)
nodes2 = getNodes(fem);
nodes2 = nodes2 + U;

[P, ~]= getInterpolationMatrix(fem, landmark_US);
results(1).landmark_US_transformed = P*nodes2;

error = results(1).landmark_US_transformed - landmark_MR;
error = sum(error.*error,2);
error = sqrt(error);
mean_error = mean(error);
sd_error = sqrt(var(error));
results(1).U = U;
results(1).error = error;

%% CPD + FEM
lambda = 0.1;
beta2_2 = 3.5;

tic;
[TY_2, V_2, cpd_results(1).sigma2, fem_2, U_2, Phi_2] = cpd_plus_boundary_fem(X, Y, Yf, lambda, beta2_2, w, errtol, maxiters, sigma2, E, nu );
time=toc; fprintf('Registration time is %f seconds\n', time);
fprintf('Registration sigma2 is %f\n', cpd_results(1).sigma2);

%% Landmark error

% original FEM
setSurfaceMesh(fem_2, Yf);

% new FEM with moved node locations (for reverse interpolation)
nodes2 = getNodes(fem_2);
elems2 = getElements(fem_2);
nodes2 = nodes2 + U_2;

[P_2, in]= getInterpolationMatrix(fem_2, landmark_US);
cpd_results(1).landmark_US_transformed = P_2*nodes2;

error_2 = cpd_results(1).landmark_US_transformed - landmark_MR;
error_2 = sum(error_2.*error_2,2);
error_2 = sqrt(error_2);
cpd_results(1).error = error_2;
cpd_results(1).U = U_2;

%% Remove data

rd = [0.025, 0.05, 0.075, 0.1 0.15 0.2 0.25];

for i=1:length(rd)
    % cut data
    XX = cut_data(X, rd(i), rd(i), [0 0 1]);
    [TY, V, results(i+1).sigma2, ~, fem, U, Phi] = cpd_fem_norigid(XX, Y, Yf, w, errtol, maxiters, sigma2, beta, E, nu);
    setSurfaceMesh(fem, Yf);
    
    nodes = getNodes(fem);
    nodes = nodes + U;
    [P, ~]= getInterpolationMatrix(fem, landmark_US);
    results(i+1).landmark_US_transformed = P*nodes;

    error = results(i+1).landmark_US_transformed - landmark_MR;
    error = sum(error.*error,2);
    error = sqrt(error);

    results(i+1).U = U;
    results(i+1).error = error;
    
    
    [~, ~, cpd_results(i+1).sigma2, fem, U, ~] = cpd_plus_boundary_fem(XX, Y, Yf, lambda, beta2_2, w, errtol, maxiters, sigma2, E, nu );
    setSurfaceMesh(fem, Yf);
    
    nodes = getNodes(fem);
    nodes = nodes + U;
    [P, ~]= getInterpolationMatrix(fem, landmark_US);
    cpd_results(i+1).landmark_US_transformed = P*nodes;

    error = cpd_results(i+1).landmark_US_transformed - landmark_MR;
    error = sum(error.*error,2);
    error = sqrt(error);

    cpd_results(i+1).U = U;
    cpd_results(i+1).error = error;
    
    fprintf('Done doing %d, errors: %f %f\n', i, mean(results(i+1).error), mean(cpd_results(i+1).error));
end

save([testcase,'.mat'], 'results', 'cpd_results');