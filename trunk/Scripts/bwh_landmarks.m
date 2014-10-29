function [fem,Phi] = bwh_landmarks(caseId,usePartialData,useRigid)
%clear; clc; close all;

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

dataPath = [root '/data/BWHTestData'];
% casePath = [ dataPath '/Case' caseId];

[landmarksMR,landmarksUS] = readBWHlandmarks(dataPath, caseId);

fprintf('Landmarks read:\n')
fprintf('MR landmarks:\n');
disp(landmarksMR);
fprintf('US landmarks:\n');
disp(landmarksUS);

fixedModelName = [ dataPath '/Case' caseId '/Input/case' caseId '-US-smooth.ply'];
fixedPartialModelName = [ dataPath '/Case' caseId '/Input/case' caseId '-US-smooth-cut10.ply'];
movingModelName = [ dataPath '/Case' caseId '/Input/case' caseId '-MR-smooth.ply'];

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

% MR Model
[nodes, ~, elems] = tetgen_mex(moving.vertices', moving.faces',[],'');
fem = fem_model(nodes', elems');
Phi = getInterpolationMatrix(fem, landmarksMR);
fprintf('MR Interpolation matrix:\n');
disp(Phi);

% TODO: rigid alignment is not robust with partial data; need to do rigid
% initialization using full data, and then FEM on partial

% display landmarks
figure(2);
subplot(1,2,1);
patch(moving,'FaceColor','blue','FaceAlpha',0.2, 'EdgeColor', 'blue', 'EdgeAlpha',0.2);
hold('on');
plot3(landmarksMR(:,1), landmarksMR(:,2), landmarksMR(:,3), 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
title('MR (moving)');
subplot(1,2,2);
patch(fixed,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
hold('on');
plot3(landmarksUS(:,1), landmarksUS(:,2), landmarksUS(:,3), 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
title('US (fixed)');

end


function [MRl, USl] = readBWHlandmarks(path,caseId)

MRfileName = [path '/Case' caseId '/Input/MR-fiducials.fcsv'];
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

USfileName = [path '/Case' caseId '/Input/US-fiducials.fcsv'];
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
