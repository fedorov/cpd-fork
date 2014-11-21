add_bcpd_paths;
% matlabpool open;

%% load data
[X.vertices,X.faces] = readPolyDataInVTK('..\data\prostate\P1.vtk');

%% create model
[nodes, ~, elems] = tetgen_mex(X.vertices', X.faces',[],'');
fem = fem_model(nodes', elems');
smesh_tree = smesh_bvtree(X.vertices', X.faces');

%% Containing element

[R, c, w] = tight_box(getNodes(fem));

% 1000 pnts in [0 1]^3
NP = 1000;
pnts = rand(NP, 3);
pnts = (pnts - 0.5).*repmat(w, [size(pnts,1) 1]);
pnts = pnts*R + repmat(c, [size(pnts,1) 1]);

% locate points
disp(['Searching for ', num2str(NP),' random points: ']);
disp('   All at once... ');
tic;
[idx, elem, N] = findContainingElement(fem, pnts);
findtime = toc;
disp(['      ...', num2str(findtime), ' (s)']);

outIdxs = find(idx == 0);
inIdxs = find(idx ~= 0);

fprintf('      Verifying...');
P = getInterpolationMatrix(fem, pnts(inIdxs,:));
nodes = getNodes(fem);
err = P*nodes - pnts(inIdxs,:);
err = sum(err.*err,2);
if (max(abs(err)) < 1e-14)
    fprintf(' PASSED :)\n');
else 
    fprintf(' FAILED :( :(\n');
end

disp('   Individually... ');
tic;
iidx = zeros(NP,1);
ielem{NP,1} = [];
iN{NP,1} = [];
for i=1:NP
    [iidxa, ielema, iNa] = findContainingElement(fem, pnts(i,:));
    iidx(i) = iidxa;
    ielem{i} = ielema{1};
    iN{i} = iNa{1};
end
findtime = toc;
disp(['      ...', num2str(findtime), ' (s)']);

%% Tree comparisons

in = is_inside(smesh_tree, pnts);
delete(smesh_tree);
