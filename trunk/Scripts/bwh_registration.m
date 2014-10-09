function [fem,u] = registerBWHcase(caseId)
clear; clc; close all;
root = '/Users/fedorov/github/cpd/trunk';

caseId = '10';

add_bcpd_paths;

dataPath = '/Users/fedorov/Documents/Projects/BRP/MR-US-registration';
casePath = [ dataPath '/Case' caseId];

[landmarksMR,landmarksUS] = readBWHlandmarks(dataPath, caseId);

fprintf('Landmarks read:')
landmarksMR
landmarksUS

fixedModelName = [ dataPath '/Case' caseId '/SmoothReg/Models/case' caseId '-US-simplified.ply'];
movingModelName = [ dataPath '/Case' caseId '/SmoothReg/Models/case' caseId '-MR-simplified.ply'];
registeredPath = [ dataPath '/Case' caseId '/SmoothReg/Models/'];

%% Read in the surfaces
[fixedVertices,fixedFaces] = read_ply(fixedModelName);
[movingVertices,movingFaces] = read_ply(movingModelName);
fprintf('Read input surfaces');

% The Slicer model has too many vertices and faces. I need to downsample
% it to use it.
numberOfFaces = 1800;

fixed.faces = fixedFaces;
fixed.vertices = fixedVertices;

moving.faces = movingFaces;
moving.vertices = movingVertices;

% fixed points
%nfx = reducepatch(pX,numberOfFaces);
% moving surface
%nfy = reducepatch(pY,numberOfFaces);

%% Non-rigid registration parameters
w=0.00; errtol=1e-4; maxiters=500; sigma2=10;
beta=0.12; E=4.8; nu=0.49;

tic;
%[TYfem, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(pX.vertices, pY.vertices, pY.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu, [], [], []);
%[TYrigidfem, ~, ~, ~, ~, ~, ~, ~] = cpd_rigid_fem(pX.vertices, [], pY.vertices, pY.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu);
%[TYrigid,~,~,~,~] = cpd_rigid(pX.vertices, pY.vertices, w, errtol, maxiters, [], [], 0, sigma2);
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
affineResultFilename = [ dataPath '/Case' caseId '/SmoothReg/Models/case' caseId '_affine.mat'];
save(affineResultFilename, 'affineResult');
affineSurfaceFilename = [ dataPath '/Case' caseId '/SmoothReg/Models/case' caseId '_affine.ply'];
write_ply(TYaffine, moving.faces, affineSurfaceFilename);

writeBWHlandmarks([ dataPath '/Case' caseId '/SmoothReg/Models/'], caseId, 'affineReg', landmarksAffine);

fprintf('Affine results saved');

landmarksAffine

tic;
[TYfem, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(fixed.vertices, TYaffine, moving.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu, [], [], []);
time=toc; fprintf('FEM registration time is %f seconds\n', time);

femResult.TY = TYfem;
femResult.sigma2 = newSigma2;
femResult.fem = fem;
femResult.u = u;
femResultFilename = [ dataPath '/Case' caseId '/SmoothReg/Models/case' caseId '_fem.mat'];
save(femResultFilename, 'femResult');
femSurfaceFilename = [ dataPath '/Case' caseId '/SmoothReg/Models/case' caseId '_fem.ply'];
write_ply(TYfem, moving.faces, femSurfaceFilename);

Phi = getInterpolationMatrix(fem, landmarksAffine);
time=toc; fprintf('FEM interpolation time is %f seconds\n', time);
fprintf('FEM results saved');

landmarksFEM = landmarksAffine+Phi*u

landmarksUS

writeBWHlandmarks([ dataPath '/Case' caseId '/SmoothReg/Models/'], caseId, 'femReg', landmarksFEM);

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
