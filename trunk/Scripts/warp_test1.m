add_bcpd_paths;

%% Read the 3D MR segmentation
[Y,fY] = readPolyDataInVTK('..\data\prostate\P069_MR_Seg.vtk');
[X,fX] = readPolyDataInVTK('..\data\prostate\P069_US_Seg.vtk');

%% Perform rigid+FEM registration
s=1; R=eye(3); t=[0;0;0];
w=0.0; errtol=1e-10; maxiters=100; sigma2=[];
beta=0.1; E=4.8; nu=0.49;

tic;
[TY, R, t, s, ~, ~, ~, fem, u, ~] = cpd_fem(X, Y, fY, w, errtol, maxiters, R, t, s, sigma2, beta, E, nu, [], [], []);
time=toc; fprintf('FEM registration time is %f seconds\n', time);

%% Write out the registered surface
writePolyDataInVTK(TY,fY,'..\data\prostate\P069_MR_Reg_Seg.vtk');

%% Read the US image which 
[US,USor,USsp] = readImageDataInVTK('..\data\prostate\P069_US_Volume.vtk');

%% Create a bounding box around the US segmentation
maxX = max(X,[],1);
minX = min(X,[],1);

offset = minX';

ds = [2; 2; 2];
ds = USsp.*ds;

% Create grid for this box
xd1 = minX(1):ds(1):maxX(1);
xd2 = minX(2):ds(2):maxX(2);
xd3 = minX(3):ds(3):maxX(3);

[xd1US,xd2US,xd3US] = meshgrid(xd1, xd2, xd3);
xdUS = [reshape(xd1US,1,[])' reshape(xd2US,1,[])' reshape(xd3US,1,[])'];

minX = USor;
maxX = USor + USsp.*(size(US)'-[1;1;1]);

x1 = minX(1):USsp(1):maxX(1);
x2 = minX(2):USsp(2):maxX(2);
x3 = minX(3):USsp(3):maxX(3);

[x1US,x2US,x3US] = meshgrid(x1, x2, x3);

order = [2 1 3];
US = permute(US,order);

tic;
dUS = interp3(x1US,x2US,x3US,US,xdUS(:,1),xdUS(:,2),xdUS(:,3));
time=toc; fprintf('Ultrasound down-sampling time is %f seconds\n', time);

offset = offset - USor;
dOr = USor;
dOr = dOr + offset;

dUS = reshape(dUS, size(xd1US));
dUS = permute(dUS,order);
writeImageDataInVTK(dUS,dOr,'..\data\prostate\P069_US_Down10.vtk',ds);
clear US; clear dUS; clear x1US; clear x2US; clear x3US;

%% Resample the ground truth registered MR
[MR,MRor,MRsp] = readImageDataInVTK('..\data\prostate\P069_MR_Reg_Volume.vtk');

minY = MRor;
maxY = MRor + MRsp.*(size(MR)'-[1;1;1]);

x1 = minY(1):MRsp(1):maxY(1);
x2 = minY(2):MRsp(2):maxY(2);
x3 = minY(3):MRsp(3):maxY(3);

[x1MR,x2MR,x3MR] = meshgrid(x1, x2, x3);
MR = permute(MR,order);

tic;
gMR = interp3(x1MR,x2MR,x3MR,MR,xdUS(:,1),xdUS(:,2),xdUS(:,3));
time=toc; fprintf('MR ground truth down-sampling time is %f seconds\n', time);

gMR = reshape(gMR, size(xd1US) );
gMR = permute(gMR,order);
writeImageDataInVTK(gMR,dOr,'..\data\prostate\P069_MR_Down10.vtk',ds);
clear MR; clear gMR; clear x1MR; clear x2MR; clear x3MR;

%% Rotate and translate the FEM using the result of rigid registration
% and get the interpolation matrix at nodes
nodes = getNodes(fem);
elems = getElements(fem);

% Write out the mesh before projection
% faces = [];

% Put elements in a matrix
% for i=1:numel(elems)
%     i/numel(elems)
%     faces = [faces; combntns(elems{i}.getNodeIdxs,3)];
%     faces = unique(faces,'rows');
% end
% writePolyDataInVTK(nodes,faces,'..\data\prostate\P069_MR_FEM.vtk');

nodes = nodes + u;
nodes = bsxfun(@plus, nodes*R'*s, t');
fem_US = fem_model( nodes, elems );
% writePolyDataInVTK(nodes,faces,'..\data\prostate\P069_US_FEM.vtk');

tic;
Phi = getInterpolationMatrix(fem, xdUS);
time=toc; fprintf('FEM interpolation time is %f seconds\n', time);

R_back = R';
s_back = 1/s;
t_back = -(R')*t/s;
u_back = -bsxfun(@plus, u*(R')*s, t');

xrUS = bsxfun(@plus, xdUS*(R_back')*s_back, t_back');

Q = Phi*u_back;
xdUS = Q + xdUS;
xdUS = xdUS*(R_back')*s_back;
%% Read the MR image
[MR,MRor,MRsp] = readImageDataInVTK('..\data\prostate\P069_MR_Volume.vtk');

%% Create the regular grid that corresponds to the MR
minY = MRor;
maxY = MRor + MRsp.*(size(MR)'-[1;1;1]);

x1 = minY(1):MRsp(1):maxY(1);
x2 = minY(2):MRsp(2):maxY(2);
x3 = minY(3):MRsp(3):maxY(3);

[x1MR,x2MR,x3MR] = meshgrid(x1, x2, x3);
MR = permute(MR,order);

%% Interpolate xUS grid on the MR image
tic;
vUS = interp3(x1MR,x2MR,x3MR,MR,xdUS(:,1),xdUS(:,2),xdUS(:,3));
time=toc; fprintf('Intensity interpolation time is %f seconds\n', time);

vUS = reshape(vUS, size(xd1US) );
vUS = permute(vUS,order);
writeImageDataInVTK(vUS,dOr,'..\data\prostate\P069_MR_DefReg10.vtk',ds);

tic;
vUS = interp3(x1MR,x2MR,x3MR,MR,xrUS(:,1),xrUS(:,2),xrUS(:,3));
time=toc; fprintf('Intensity interpolation time is %f seconds\n', time);

vUS = reshape(vUS, size(xd1US) );
vUS = permute(vUS,order);
writeImageDataInVTK(vUS,dOr,'..\data\prostate\P069_MR_RigidReg10.vtk',ds);