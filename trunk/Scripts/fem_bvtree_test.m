add_bcpd_paths;
% matlabpool open;

%% beam model construction

beam = fem_beam_model([2 2 2], [16, 16, 16]);
% duplicate in a regular fem_model
model = fem_model( getNodes(beam), getElements(beam) );

%% Material parameters
E = 480; nu = 1/3;
D = fem_material_linear(E, nu);

%% Containing element

[R, c, w] = tight_box(getNodes(model));

% 1000 pnts in [0 1]^3
NP = 100000;
pnts = rand(NP, 3);
pnts = (pnts - 0.5).*repmat(w, [size(pnts,1) 1]);
pnts = pnts*R + repmat(c, [size(pnts,1) 1]);

% locate points
%[eidx, elem, eN] = findContainingElement(beam, pnts);
disp(['Searching for ', num2str(NP),' random points: ']);
disp('   All at once... ');
tic;
[idx, elem, N] = findContainingElement(model, pnts);
findtime = toc;
disp(['      ...', num2str(findtime), ' (s)']);
fprintf('      Verifying...');
P = getInterpolationMatrix(fem, pnts);
nodes = getNodes(fem);
err = P*nodes - pnts;
err = sum(err.*err,2);
if (max(abs(err)) > 1e-14)
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
    [iidxa, ielema, iNa] = findContainingElement(model, pnts(i,:));
    iidx(i) = iidxa;
    ielem{i} = ielema{1};
    iN{i} = iNa{1};
end
findtime = toc;
disp(['      ...', num2str(findtime), ' (s)']);

% 
% disp('   In parallel ... ');
% tic;
% pidx = zeros(NP,1);
% pelem{NP,1} = [];
% pN{NP,1} = [];
% for i=1:NP
%     [pidxa, pelema, pNa] = findContainingElement(model, pnts(i,:));
%     pidx(i) = pidxa;
%     pelem{i} = pelema{1};
%     pN{i} = pNa{1};
% end
% findtime = toc;
% disp(['      ...', num2str(findtime), ' (s)']);
% 
% matlabpool close;