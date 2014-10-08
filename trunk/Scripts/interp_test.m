add_bcpd_paths;

%% Read the 3D segmentations
%[MR.X,MR.F] = readPolyDataInVTK('..\data\prostate\P069_MR_Seg.vtk');
%[US.X,US.F] = readPolyDataInVTK('..\data\prostate\P069_US_Seg.vtk');
[X, F] = read_obj('../data/synthetic/oct2.obj');

%% Generate tet meshes
[nodes, ~, elems] = tetgen_mex(X', F',[],'');
fem = fem_model(nodes', elems');
setAbortJacobian(fem, 1e-15);      % abort on inverted elements
setSurfaceMesh(fem, F);

%% Interpolate
Npoints = 10000;
[ R, c, w, B ] = tight_box( X );

pnts = rand(Npoints, 3);
pnts = (pnts - 0.5).*repmat(1.25*w, [size(pnts,1) 1]);
pnts = pnts*R + repmat(c, [size(pnts,1) 1]);

tic;
[idx1,~,~] = findContainingElement(fem, pnts);
time=toc; fprintf('Find Containing element fast is %f seconds\n', time);
fprintf('Verifying...');
P = getInterpolationMatrix(fem, pnts);
nodes = getNodes(fem);
err = P*nodes - pnts;
err = sum(err.*err,2);
if (max(abs(err)) < 1e-14)
    fprintf(' PASSED :)\n');
else 
    fprintf(' FAILED :( :(\n');
end

tic;
[idx2,~,~] = findContainingElementSlow(fem, pnts);
time=toc; fprintf('Find Containing element slow is %f seconds\n', time);

if ( sum(abs(idx1-idx2)) > 0 ) 
    display('findContainingElement does NOT return consistent results.');
    disp('Error at locations: ');
    disp(find(abs(idx1-idx2) > 0));
else 
    display('findContainingElement consistent :) :).');
end

