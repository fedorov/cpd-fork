add_bcpd_paths;

%% beam model construction

beam = fem_beam_model([6 8 10], [10, 12, 16]);

% duplicate in a regular fem_model
model = fem_model( getNodes(beam), getElements(beam) );

%% Material parameters

E = 480; nu = 1/3;
D = fem_material_linear(E, nu);

%% Stiffness matrix computations

fprintf('Computing stiffness matrix from model\n');
tic;
K = getStiffnessMatrix(beam, D, 1);
model_time = toc;
fprintf('Elapsed time is %f seconds\n', model_time);
fprintf('...again\n');
tic;
K2 = getStiffnessMatrix(beam, D, 1);
model_time_2 = toc;
fprintf('Elapsed time is %f seconds\n', model_time_2);

fprintf('Computing stiffness matrix raw\n');
tic;
K3 = fem_stiffness_matrix(getNodes(beam),getElements(beam),D);
model_time_3 = toc;
fprintf('Elapsed time is %f seconds\n', model_time_3);

%% Interpolation matrix

nodes = getNodes(beam);
fprintf('Computing interpolation matrix\n');
tic;
P = getInterpolationMatrix(beam, nodes);
interpolation_time = toc;
fprintf('Elapsed time is %f seconds\n', interpolation_time);

%% Containing element

[eidx, elem, eN] = findContainingElement(beam, [0 0 0]);