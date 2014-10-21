function [fem,u] = bwh_registration(caseId,usePartialData)
%clear; clc; close all;
root = '/Users/fedorov/github/cpd/trunk';

usePartialData = 1;
caseId='9';

add_bcpd_paths;

dataPath = '/Users/fedorov/Documents/Projects/BRP/MR-US-registration';
casePath = [ dataPath '/Case' caseId];

[landmarksMR,landmarksUS] = readBWHlandmarks(dataPath, caseId);

fprintf('Landmarks read:')
landmarksMR
landmarksUS

fixedModelName = [ dataPath '/Case' caseId '/SmoothReg/case' caseId '-US-smooth.ply'];
fixedPartialModelName = [ dataPath '/Case' caseId '/SmoothReg/case' caseId '-US-smooth-cut10.ply'];
movingModelName = [ dataPath '/Case' caseId '/SmoothReg/case' caseId '-MR-smooth.ply'];
registeredPath = [ dataPath '/Case' caseId '/CPD_registration/'];

%% Read in the surfaces
[fixedVertices,fixedFaces] = read_ply(fixedModelName);
[fixedPartialVertices,fixedPartialFaces] = read_ply(fixedPartialModelName);
[movingVertices,movingFaces] = read_ply(movingModelName);
fprintf('Read input surfaces');

% The Slicer model has too many vertices and faces. I need to downsample
% it to use it.
numberOfFaces = 1800;

fixed.faces = fixedFaces;
fixed.vertices = fixedVertices;

fixedPartial.faces = fixedPartialFaces;
fixedPartial.vertices = fixedPartialVertices;

moving.faces = movingFaces;
moving.vertices = movingVertices;

% fixed points
%nfx = reducepatch(pX,numberOfFaces);
% moving surface
%nfy = reducepatch(pY,numberOfFaces);

%% Non-rigid registration parameters
w=0.00; errtol=1e-4; maxiters=500; sigma2=10;
beta=0.12; E=4.8; nu=0.49;

% TODO: rigid alignment is not robust with partial data; need to do rigid
% initialization using full data, and then FEM on partial

tic;
%[TYfem, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(pX.vertices, pY.vertices, pY.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu, [], [], []);
%[TYrigidfem, ~, ~, ~, ~, ~, ~, ~] = cpd_rigid_fem(pX.vertices, [], pY.vertices, pY.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu);

%[TYrigid,B_affine,t_affine,~,~] = cpd_rigid(fixed.vertices, moving.vertices, w, errtol, maxiters, [], [], 0, sigma2);

[TYaffine,B_affine,t_affine,~,~] = cpd_affine(fixed.vertices, moving.vertices, w, errtol, maxiters, [], [], sigma2);
time=toc; fprintf('Affine registration time is %f seconds\n', time);

landmarksAffine = [];
numLandmarks = size(landmarksMR,1);
for l=1:numLandmarks
    lm = landmarksMR(l,:);
    registered = [];
    registered = bsxfun(@plus,lm*(B_affine'),t_affine');
    landmarksAffine = [landmarksAffine; registered];
end

landmarksAffine

affineResult.TY = TYaffine;
affineResult.B = B_affine;
affineResult.t = t_affine;
affineResultFilename = [ registeredPath '/case' caseId '_affine.mat'];
save(affineResultFilename, 'affineResult');
affineSurfaceFilename = [ registeredPath '/case' caseId '_affine.ply'];
write_ply(TYaffine, moving.faces, affineSurfaceFilename);

writeAffineTransform([ registeredPath '/case' caseId '_affine.tfm'], affineResult.B, affineResult.t);

writeBWHlandmarks(registeredPath, caseId, 'affineReg', landmarksAffine);

fprintf('Affine results saved');

landmarksAffine

% estimate of the missing data
w=0.20;

tic;
femType = '';
if usePartialData
    femType = 'Partial';
    fprintf('Registering with partial data:')
    size(fixedPartial.vertices)
    [TYfem, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(fixedPartial.vertices, TYaffine, moving.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu, [], [], []);
else
    femType = 'Full';
    [TYfem, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(fixed.vertices, TYaffine, moving.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu, [], [], []);
end


time=toc; fprintf('FEM registration time is %f seconds\n', time);

femResult.TY = TYfem;
femResult.sigma2 = newSigma2;
femResult.fem = fem;
femResult.u = u;
femResultFilename = [ registeredPath '/case' caseId '_fem_' femType '.mat'];
save(femResultFilename, 'femResult');
femSurfaceFilename = [ registeredPath '/case' caseId '_fem_' femType '.ply'];
write_ply(TYfem, moving.faces, femSurfaceFilename);
meshFileName = [ registeredPath '/case' caseId '_fem_result_mesh_' femType '.vtk']
writeFEMvtk(fem,u,meshFileName,1)

Phi = getInterpolationMatrix(fem, landmarksAffine);
time=toc; fprintf('FEM interpolation time is %f seconds\n', time);
fprintf('FEM results saved');

landmarksFEM = landmarksAffine+Phi*u

landmarksUS

writeBWHlandmarks(registeredPath, caseId, ['femReg_' femType], landmarksFEM);

initial_error = sqrt(sum((landmarksMR-landmarksUS).*(landmarksMR-landmarksUS),2));
affine_error = sqrt(sum((landmarksAffine-landmarksUS).*(landmarksAffine-landmarksUS),2));
fem_error = sqrt(sum((landmarksFEM-landmarksUS).*(landmarksFEM-landmarksUS),2));

end


function [MRl, USl] = readBWHlandmarks(path,caseId)

MRfileName = [path '/Annotation/Case' caseId '/MR-fiducials.fcsv'];
C=textscan(fopen(MRfileName,'r'),'%s','Delimiter','\n');
num=size(C{1},1);
MRl = [];
for i=4:num
  coordStr = strsplit(C{1}{i},',');
  coordX = coordStr(2);
  coordY = coordStr(3);
  coordZ = coordStr(4);
  MRl = [ MRl; str2num(coordX{1}) str2num(coordY{1}) str2num(coordZ{1}) ];
end

USfileName = [path '/Annotation/Case' caseId '/US-fiducials.fcsv'];
C=textscan(fopen(USfileName,'r'),'%s','Delimiter','\n');
num=size(C{1},1);
USl = [];
for i=4:num
  coordStr = strsplit(C{1}{i},',');
  coordX = coordStr(2);
  coordY = coordStr(3);
  coordZ = coordStr(4);
  USl = [ USl; str2num(coordX{1}) str2num(coordY{1}) str2num(coordZ{1}) ];
end

end

function writeBWHlandmarks(path,caseId,name,L)

fileName = [path '/case' caseId '-' name '.fcsv'];
fprintf('Writing landmarks to %s\n',fileName);
fid = fopen(fileName,'w');
fprintf(fid,'# Markups fiducial file version = 4.3');
fprintf(fid,'# CoordinateSystem = 0');
fprintf(fid,'# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID');

for i=1:size(L,1)
  fprintf(fid,'vtkMRMLMarkupsFiducialNode_%i,%f,%f,%f,0,0,0,1,1,1,0,%s%i,,',i,L(i,1),L(i,2),L(i,3),name,i);
end

fclose(fid);

end

function writeAffineTransform(filename, B, t)
tMatrix = zeros(4);
tMatrix(1:3,1:3) = B;
tMatrix(:,4) = [t',1];
lps2ras = eye(4);
lps2ras(1,1) = -1;
lps2ras(2,2) = -1;
ras2lps = lps2ras;
tMatrix = inv(lps2ras*tMatrix*ras2lps);


fid=fopen(filename, 'w');

fprintf(fid,'#Insight Transform File V1.0\n');
fprintf(fid,'#Transform 0\n');
fprintf(fid,'Transform: AffineTransform_double_3_3\n');
fprintf(fid,'Parameters: ');

for i=1:3
    for j=1:3
        fprintf(fid,'%f ',tMatrix(i,j));
    end
end

for j=1:3
    fprintf(fid,'%f ',tMatrix(j,4));
end

fprintf(fid,'\n');
fprintf(fid,'FixedParameters: 0 0 0\n');

fclose(fid);

end