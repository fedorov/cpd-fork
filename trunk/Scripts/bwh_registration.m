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

fprintf('FEM Prefix: ');
disp(femPrefix);

% caseId='9';
caseId=num2str(caseId);

add_bcpd_paths;

dataPath = '/Users/fedorov/Documents/Projects/BRP/MR-US-registration';
% casePath = [ dataPath '/Case' caseId];
% dataPath = [root '/data/BWHTestData'];
% casePath = [ dataPath '/Case' caseId];

[landmarksMR,landmarksUS] = readBWHlandmarks(dataPath, caseId);

fprintf('Landmarks read:\n')
fprintf('MR landmarks:\n');
disp(landmarksMR);
fprintf('US landmarks:\n');
disp(landmarksUS);

fixedModelName = [ dataPath '/Case' caseId '/SmoothReg/case' caseId '-US-smooth.ply'];
fixedPartialModelName = [ dataPath '/Case' caseId '/SmoothReg/case' caseId '-US-smooth-cut10.ply'];
movingModelName = [ dataPath '/Case' caseId '/SmoothReg/case' caseId '-MR-smooth.ply'];
registeredPath = [ dataPath '/Case' caseId '/CPD_registration/'];
% fixedModelName = [ dataPath '/Case' caseId '/Input/case' caseId '-US-smooth.ply'];
% fixedPartialModelName = [ dataPath '/Case' caseId '/Input/case' caseId '-US-smooth-cut10.ply'];
% movingModelName = [ dataPath '/Case' caseId '/Input/case' caseId '-MR-smooth.ply'];
% registeredPath = [ dataPath '/Case' caseId '/Results2/'];

%% Read in the surfaces
[fixedVertices,fixedFaces] = read_ply(fixedModelName);
[fixedPartialVertices,fixedPartialFaces] = read_ply(fixedPartialModelName);
[movingVertices,movingFaces] = read_ply(movingModelName);
fprintf('Read input surfaces\n');

fixed.faces = fixedFaces;
fixed.vertices = fixedVertices;

fixedPartial.faces = fixedPartialFaces;
fixedPartial.vertices = fixedPartialVertices;

moving.faces = movingFaces;
moving.vertices = movingVertices;

% The Slicer model has too many vertices and faces. I need to downsample
% it to use it.
% numberOfFaces = 1800;
% fixed points
%nfx = reducepatch(pX,numberOfFaces);
% moving surface
%nfy = reducepatch(pY,numberOfFaces);

%% Non-rigid registration parameters
errtol=1e-4; maxiters=500; sigma2=10;
beta=0.3; E=4.8; nu=0.49;
if usePartialData
    w = 0.1;  % 10% missing data?
else
    w = 0.01;  % maybe 1% outliers
end

% TODO: rigid alignment is not robust with partial data; need to do rigid
% initialization using full data, and then FEM on partial

tic;
%[TYfem, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(pX.vertices, pY.vertices, pY.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu, [], [], []);
%[TYrigidfem, ~, ~, ~, ~, ~, ~, ~] = cpd_rigid_fem(pX.vertices, [], pY.vertices, pY.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu);

if useRigid == 1
    fprintf('Before rigid registration\n');
    [TYrigidOrAffine,B_rigid,t_rigid,s_rigid,~] = cpd_rigid(fixed.vertices, moving.vertices, w, errtol, maxiters);
    B_rigid = B_rigid*s_rigid;  % fold in scale
    fprintf('Rigid registration completed\n');
    landmarksRigid = bsxfun(@plus,landmarksMR*(B_rigid'),t_rigid');

    rigidResult = struct('TY', TYrigidOrAffine, 'B', B_rigid, 't', t_rigid);
    disp([rigidResult.B rigidResult.t]);

    % write results
    rigidResultFilename = [ registeredPath '/case' caseId '_rigid.mat'];
    save(rigidResultFilename, 'rigidResult');
    rigidSurfaceFilename = [ registeredPath '/case' caseId '_rigid.ply'];
    write_ply(TYrigidOrAffine, moving.faces, rigidSurfaceFilename);
    writeLinearTransform([ registeredPath '/case' caseId '_rigid.tfm'], rigidResult.B, rigidResult.t);
    writeBWHlandmarks(registeredPath, caseId, 'rigidReg', landmarksRigid);

    fprintf('Rigid results saved\n');

    landmarksPreRegistered = landmarksRigid;
else
 [TYrigidOrAffine,B_affine,t_affine,~,~] = cpd_affine(fixed.vertices, moving.vertices, w, errtol, maxiters);
    landmarksAffine = bsxfun(@plus,landmarksMR*(B_affine'),t_affine');
    affineResult = struct('TY', TYrigidOrAffine, 'B', B_affine,'t', t_affine);
    fprintf('Affine results saved\n');
    disp([affineResult.B affineResult.t]);

    % write results
    affineResultFilename = [ registeredPath '/case' caseId '_affine.mat'];
    save(affineResultFilename, 'affineResult');
    affineSurfaceFilename = [ registeredPath '/case' caseId '_affine.ply'];
    write_ply(TYrigidOrAffine, moving.faces, affineSurfaceFilename);
    writeLinearTransform([ registeredPath '/case' caseId '_affine.tfm'], affineResult.B, affineResult.t);
    writeBWHlandmarks(registeredPath, caseId, 'affineReg', landmarksAffine);

    fprintf('Affine results saved\n');
    
    landmarksPreRegistered = landmarksAffine;

end

time=toc; fprintf('Initial registration time is %f seconds\n', time);

tic;
if usePartialData
    fprintf('Registering with partial data:\n')
    size(fixedPartial.vertices)
    [TYfem, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(fixedPartial.vertices, TYrigidOrAffine, moving.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, [], beta, E, nu);
else
    [TYfem, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(fixed.vertices, TYrigidOrAffine, moving.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, [], beta, E, nu);
end


time=toc; fprintf('FEM registration time is %f seconds\n', time);

femResult.TY = TYfem;
femResult.sigma2 = newSigma2;
femResult.fem = fem;
femResult.u = u;

% write results
femResultFilename = [ registeredPath '/case' caseId '_' femPrefix '_fem.mat'];
save(femResultFilename, 'femResult');
femSurfaceFilename = [ registeredPath '/case' caseId '_' femPrefix '_fem.ply'];
write_ply(TYfem, moving.faces, femSurfaceFilename);
meshFileName = [ registeredPath '/case' caseId '_' femPrefix '_fem_result_mesh.vtk'];
writeFEMvtk(fem,u,meshFileName,1);
time=toc; fprintf('FEM interpolation time is %f seconds\n', time);
fprintf('FEM results saved\n');

% Interpolate landmarks
Phi = getInterpolationMatrix(fem, landmarksPreRegistered);
landmarksFEM = landmarksPreRegistered+Phi*u;

fprintf('FEM landmarks:\n');
disp(landmarksFEM);
fprintf('US landmarks\n');
disp(landmarksUS);

initial_error = sqrt(sum((landmarksMR-landmarksUS).*(landmarksMR-landmarksUS),2));
affine_error = sqrt(sum((landmarksPreRegistered-landmarksUS).*(landmarksPreRegistered-landmarksUS),2));
fem_error = sqrt(sum((landmarksFEM-landmarksUS).*(landmarksFEM-landmarksUS),2));

fprintf('Initial error:\n');
disp(initial_error);
fprintf('Rigid/Affine error:\n');
disp(affine_error);
fprintf('FEM error:\n');
disp(fem_error);

% Write landmarks and errors
writeBWHlandmarks(registeredPath, caseId, ['_' femPrefix 'femReg_'], landmarksFEM);
writeBWHErrors(initial_error, affine_error, fem_error, [registeredPath '/case' caseId femPrefix '_errors.txt']);

% display landmarks
figure(2);
subplot(1,3,1);
patch(moving,'FaceColor','blue','FaceAlpha',0.2, 'EdgeColor', 'blue', 'EdgeAlpha',0.2);
hold('on');
plot3(landmarksMR(:,1), landmarksMR(:,2), landmarksMR(:,3), 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
title('MR (moving)');
subplot(1,3,2);
patch(fixed,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
hold('on');
plot3(landmarksUS(:,1), landmarksUS(:,2), landmarksUS(:,3), 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
title('US (fixed)');
subplot(1,3,3);
plot3(landmarksUS(:,1), landmarksUS(:,2), landmarksUS(:,3), 'or', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
hold('on');
plot3(landmarksFEM(:,1), landmarksFEM(:,2), landmarksFEM(:,3), 'ob', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
plot3([landmarksFEM(:,1)'; landmarksUS(:,1)'],...
    [landmarksFEM(:,2)'; landmarksUS(:,2)'],...
    [landmarksFEM(:,3)'; landmarksUS(:,3)'], 'k-');
patch(fixed,'FaceColor','red','FaceAlpha',0.05, 'EdgeColor', 'red', 'EdgeAlpha',0.05);
moved.vertices = femResult.TY;
moved.faces = moving.faces;
patch(moved,'FaceColor','blue','FaceAlpha',0.05, 'EdgeColor', 'blue', 'EdgeAlpha',0.05);
title('Landmark Error');

end


function [MRl, USl] = readBWHlandmarks(path,caseId)

MRfileName = [path '/Annotation/Case' caseId '/MR-fiducials.fcsv'];
% MRfileName = [path '/Case' caseId '/Input/MR-fiducials.fcsv'];
fprintf('Reading MR fiducials from %s\n',MRfileName);
C=textscan(fopen(MRfileName,'r'),'%s','Delimiter','\n');
num=size(C{1},1);
MRl = zeros(num-3,3);
for i=4:num
  coordStr = strsplit(C{1}{i},',');
  coordX = coordStr(2);
  coordY = coordStr(3);
  coordZ = coordStr(4);
  MRl(i-3,:) = [ str2double(coordX{1}) str2double(coordY{1}) str2double(coordZ{1}) ];
end

USfileName = [path '/Annotation/Case' caseId '/US-fiducials.fcsv'];
% USfileName = [path '/Case' caseId '/Input/US-fiducials.fcsv'];
fprintf('Reading US fiducials from %s\n',USfileName);
C=textscan(fopen(USfileName,'r'),'%s','Delimiter','\n');
num=size(C{1},1);
USl = zeros(num-3,3);
for i=4:num
  coordStr = strsplit(C{1}{i},',');
  coordX = coordStr(2);
  coordY = coordStr(3);
  coordZ = coordStr(4);
  USl(i-3,:) = [str2double(coordX{1}) str2double(coordY{1}) str2double(coordZ{1}) ];
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
