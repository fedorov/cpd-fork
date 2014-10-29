function [fem, Phi] = bwh_rigid_registration(caseId,usePartialData,useRigid)
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

% fixed points
%nfx = reducepatch(pX,numberOfFaces);
% moving surface
%nfy = reducepatch(pY,numberOfFaces);

%% Non-rigid registration parameters
w=0.05; errtol=1e-4; maxiters=500;

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

    landmarksPreRegistered = landmarksRigid;
else
    [TYrigidOrAffine,B_affine,t_affine,~,~] = cpd_affine(fixed.vertices, moving.vertices, w, errtol, maxiters);
    landmarksAffine = bsxfun(@plus,landmarksMR*(B_affine'),t_affine');
    affineResult = struct('TY', TYrigidOrAffine, 'B', B_affine,'t', t_affine);
    fprintf('Affine results saved\n');
    disp([affineResult.B affineResult.t]);
    landmarksPreRegistered = landmarksAffine;

end
time=toc; fprintf('Initial registration time is %f seconds\n', time);
        
% MR Model
[nodes, ~, elems] = tetgen_mex(TYrigidOrAffine', moving.faces',[],'');
fem = fem_model(nodes', elems');
Phi = getInterpolationMatrix(fem, landmarksPreRegistered);
fprintf('MR Interpolation matrix:\n');
disp(Phi);

if usePartialData
    X = fixedPartial.vertices;
else
    X = fixed.vertices;
end
Y = moving.vertices;
TY = TYrigidOrAffine;

        [az, el] = view;
        clf;
        if (~exist('FV','var') || isempty(FV))
            plot3(X(:,1), X(:,2), X(:,3),'.r','MarkerSize',10);
        else 
            patch(FV,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
        end
        hold on;
        % axis([-3,3,-3,3,-3,3]);
        % without deformation
        plot3(Y(:,1), Y(:,2), Y(:,3),'ob','MarkerSize', 5, 'MarkerFaceColor', 'g');
        plot3(TY(:,1), TY(:,2), TY(:,3),'ob','MarkerSize', 5, 'MarkerFaceColor', 'b');
        % axis([-3,3,-3,3,-3,3]);
        legend({'target','orig','rest'},'location','NEo')
        view(az, el);
        drawnow;
        
        
figure(2);
subplot(1,3,1);
patch(moving,'FaceColor','blue','FaceAlpha',0.2, 'EdgeColor', 'blue', 'EdgeAlpha',0.2);
hold('on');
plot3(landmarksMR(:,1), landmarksMR(:,2), landmarksMR(:,3), 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
title('MR (moving)');
subplot(1,3,2);
moved.vertices = TY;
moved.faces = moving.faces;
patch(moved,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
hold('on');
plot3(landmarksPreRegistered(:,1), landmarksPreRegistered(:,2), landmarksPreRegistered(:,3), 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
title('MR (moved)');
subplot(1,3,3);
plot3(landmarksPreRegistered(:,1), landmarksPreRegistered(:,2), landmarksPreRegistered(:,3), 'or', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
hold('on');
plot3(landmarksMR(:,1), landmarksMR(:,2), landmarksMR(:,3), 'ob', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
plot3([landmarksPreRegistered(:,1)'; landmarksMR(:,1)'],...
    [landmarksPreRegistered(:,2)'; landmarksMR(:,2)'],...
    [landmarksPreRegistered(:,3)'; landmarksMR(:,3)'], 'k-');
patch(moved,'FaceColor','red','FaceAlpha',0.05, 'EdgeColor', 'red', 'EdgeAlpha',0.05);
patch(moving,'FaceColor','blue','FaceAlpha',0.05, 'EdgeColor', 'blue', 'EdgeAlpha',0.05);
title('Landmark Error');


initial_error = sqrt(sum((landmarksMR-landmarksUS).*(landmarksMR-landmarksUS),2));
affine_error = sqrt(sum((landmarksPreRegistered-landmarksUS).*(landmarksPreRegistered-landmarksUS),2));
disp(initial_error);
disp(affine_error);

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
