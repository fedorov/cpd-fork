function [fem,u] = bwh_registration(caseId,usePartialData,useRigid)
%clear; clc; close all;
root = '/Users/fedorov/github/cpd/trunk';

femPrefix = '';
if usePartialData
    femPrefix = [femPrefix '_partial'];
else
    femPrefix = [femPrefix '_full'];
end
    
if useRigid
    femPrefix = [femPrefix '_rigid'];
else
    femPrefix = [femPrefix '_affine'];
end

femPrefix

% caseId='9';
caseId=num2str(caseId);

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


if useRigid == 1
    fprintf('Before rigid registration\n');
    [TYrigidOrAffine,B_rigid,t_rigid,~,~] = cpd_rigid(fixed.vertices, moving.vertices, w, errtol, maxiters, [], [], 0, sigma2);
    fprintf('Rigid registration completed\n');
    landmarksRigid = [];
    numLandmarks = size(landmarksMR,1);
    for l=1:numLandmarks
        lm = landmarksMR(l,:);
        registered = [];
        registered = bsxfun(@plus,lm*(B_rigid'),t_rigid');
        landmarksRigid = [landmarksRigid; registered];
    end

    rigidResult.TY = TYrigidOrAffine;
    rigidResult.B = B_rigid;
    rigidResult.t = t_rigid;
    rigidResultFilename = [ registeredPath '/case' caseId '_rigid.mat'];
    save(rigidResultFilename, 'rigidResult');
    rigidSurfaceFilename = [ registeredPath '/case' caseId '_rigid.ply'];
    write_ply(TYrigidOrAffine, moving.faces, rigidSurfaceFilename);

    writeLinearTransform([ registeredPath '/case' caseId '_rigid.tfm'], rigidResult.B, rigidResult.t);

    writeBWHlandmarks(registeredPath, caseId, 'rigidReg', landmarksRigid);

    fprintf('Rigid results saved');

    landmarksPreRegistered = landmarksRigid;
else
    [TYrigidOrAffine,B_affine,t_affine,~,~] = cpd_affine(fixed.vertices, moving.vertices, w, errtol, maxiters, [], [], sigma2);
    landmarksAffine = [];
    numLandmarks = size(landmarksMR,1);
    for l=1:numLandmarks
        lm = landmarksMR(l,:);
        registered = [];
        registered = bsxfun(@plus,lm*(B_affine'),t_affine');
        landmarksAffine = [landmarksAffine; registered];
    end

    affineResult.TY = TYrigidOrAffine;
    affineResult.B = B_affine;
    affineResult.t = t_affine;
    affineResultFilename = [ registeredPath '/case' caseId '_affine.mat'];
    save(affineResultFilename, 'affineResult');
    affineSurfaceFilename = [ registeredPath '/case' caseId '_affine.ply'];
    write_ply(TYrigidOrAffine, moving.faces, affineSurfaceFilename);

    writeLinearTransform([ registeredPath '/case' caseId '_affine.tfm'], affineResult.B, affineResult.t);

    writeBWHlandmarks(registeredPath, caseId, 'affineReg', landmarksAffine);

    fprintf('Affine results saved');
    
    landmarksPreRegistered = landmarksAffine;

end

time=toc; fprintf('Initial registration time is %f seconds\n', time);

% estimate of the missing data
w=0.2;

tic;
if usePartialData
    fprintf('Registering with partial data:')
    size(fixedPartial.vertices)
    [TYfem, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(fixedPartial.vertices, TYrigidOrAffine, moving.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu, [], [], []);
else
    [TYfem, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(fixed.vertices, TYrigidOrAffine, moving.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu, [], [], []);
end


time=toc; fprintf('FEM registration time is %f seconds\n', time);

femResult.TY = TYfem;
femResult.sigma2 = newSigma2;
femResult.fem = fem;
femResult.u = u;

femResultFilename = [ registeredPath '/case' caseId '_' femPrefix '_fem.mat'];
save(femResultFilename, 'femResult');
femSurfaceFilename = [ registeredPath '/case' caseId '_' femPrefix '_fem.ply'];
write_ply(TYfem, moving.faces, femSurfaceFilename);
meshFileName = [ registeredPath '/case' caseId '_' femPrefix '_fem_result_mesh.vtk'];
writeFEMvtk(fem,u,meshFileName,1);

Phi = getInterpolationMatrix(fem, landmarksPreRegistered);
time=toc; fprintf('FEM interpolation time is %f seconds\n', time);
fprintf('FEM results saved');

landmarksFEM = landmarksPreRegistered+Phi*u

landmarksUS

writeBWHlandmarks(registeredPath, caseId, ['_' femPrefix 'femReg_'], landmarksFEM);

initial_error = sqrt(sum((landmarksMR-landmarksUS).*(landmarksMR-landmarksUS),2))
affine_error = sqrt(sum((landmarksPreRegistered-landmarksUS).*(landmarksPreRegistered-landmarksUS),2))
fem_error = sqrt(sum((landmarksFEM-landmarksUS).*(landmarksFEM-landmarksUS),2))

writeBWHErrors(initial_error, affine_error, fem_error, [registeredPath '/case' caseId femPrefix '_errors.txt']);

end


function [MRl, USl] = readBWHlandmarks(path,caseId)

MRfileName = [path '/Annotation/Case' caseId '/MR-fiducials.fcsv'];
fprintf('Reading MR fiducials from %s\n',MRfileName);
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
fprintf('Reading US fiducials from %s\n',USfileName);
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
fprintf(fid,'# Markups fiducial file version = 4.3\n');
fprintf(fid,'# CoordinateSystem = 0\n');
fprintf(fid,'# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n');

for i=1:size(L,1)
  fprintf(fid,'vtkMRMLMarkupsFiducialNode_%i,%f,%f,%f,0,0,0,1,1,1,0,%s%i,,\n',i,L(i,1),L(i,2),L(i,3),name,i);
end

fclose(fid);

end

function writeLinearTransform(filename, B, t)
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

function writeBWHErrors(initial_error, affine_error, fem_error, filename)
fid = fopen(filename,'w');

fprintf(fid,'InitialError;');
for i=1:size(initial_error)
    fprintf(fid,'%f;',initial_error(i));
end
fprintf(fid,'\n');

fprintf(fid,'LinearRegError;');
for i=1:size(affine_error)
    fprintf(fid,'%f;',affine_error(i));
end
fprintf(fid,'\n');

fprintf(fid,'FinalRegError;');
for i=1:size(fem_error)
    fprintf(fid,'%f;',fem_error(i));
end
fprintf(fid,'\n');

end