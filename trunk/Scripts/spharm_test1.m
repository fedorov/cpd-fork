add_bcpd_paths;
clear; close all; clc;

% Load a prostate
[Y,fY] = read_ply('C:\data\ProstateQueens\P034-data\MR_model.ply');

numberOfFaces = 2000;

pY.faces = fY;
pY.vertices = Y;

nfy = reducepatch(pY,numberOfFaces);

%patch(nfy,'FaceColor','blue','FaceAlpha',0.2, 'EdgeColor', 'blue', 'EdgeAlpha',0.2);

vertices = nfy.vertices;
faces = nfy.faces;
vertices = bsxfun(@minus, vertices, mean(vertices));

save('..\data\prostate\prostate_034_obj.mat', 'faces', 'vertices' );

% SPHARM parameterization
paramCALD = struct;
paramCALD.MeshGridSize = 50; paramCALD.MaxSPHARMDegree = 6; paramCALD.Tolerance = 2;
paramCALD.Smoothing = 2; paramCALD.Iteration = 100; paramCALD.LocalIteration = 10;
paramCALD.t_major = 'z'; paramCALD.SelectDiagonal = 'ShortDiag'; paramCALD.OutDirectory = '..\data\prostate\';

'ParamCALD';
Objs = {'..\data\prostate\prostate_034_obj.mat'};

outName_parameterization = SpharmMatParameterization(paramCALD, Objs, 'ParamCALD');

% SPHARM expansion
outNames_expansion = SpharmMatExpansion(paramCALD, outName_parameterization, 'ExpLSF');

% Load the SPHARM expansion
load( '..\data\prostate\prostate_034_CALD_LSF_des.mat' );

step_rot = pi/100;
metric = [];
full_fvec = fvec;

for alpha=-pi/2:step_rot:pi/2
    rvec = spharm_rotate(alpha, 0.0, 0.0, full_fvec, paramCALD.MaxSPHARMDegree);
    metric = [metric; norm(full_fvec - rvec)];
end

alpha=-pi/2:step_rot:pi/2;
h0 = figure;

% axes1 = axes('Parent',h0,'YGrid','on',...
%     'XTickLabel',{'-pi/2','-pi/4','0','pi/4','pi/2'},...
%     'XTick',[-1.571 -0.785 0 0.785 1.571],...
%     'XGrid','on',...
%     'FontSize',20);

% xlim(axes1,[-1.571 1.571]);
% box(axes1,'on');
% hold(axes1,'all');

% Create plot
plot(alpha,metric,'LineWidth',2);

set(gca,'YGrid','on',...
    'XTickLabel',{'-pi/2','-pi/4','0','pi/4','pi/2'},...
    'XTick',[-1.571 -0.785 0 0.785 1.571],...
    'XGrid','on',...
    'FontSize',20);

% Create title
title('No Missing Data','FontSize',20);

% Create xlabel
xlabel('Rotation about the lateral axis of the prostate in radians',...
    'FontSize',20);

% Create ylabel
ylabel('SPHARM distance metric','FontSize',20);
% box(axes1,'off');
% hold(axes1,'off');
close(h0);

vertices_full = vertices;

step_w = 0.05;
for w=0.0:step_w:0.5
    vertices = cut_data( vertices_full, w/2, w/2, [0 0 1] );
    f = delaunay(vertices(:,1),vertices(:,2),vertices(:,3));
    faces = tet2tri(f,[vertices(:,1) vertices(:,2) vertices(:,3)],1);
    
    pY.faces = faces;
    pY.vertices = vertices;

    nfy = reducepatch(pY,numberOfFaces);
    
    vertices = nfy.vertices;
    faces = nfy.faces;
    
    % debug
    %figure, patch('Vertices',vertices,'Faces',faces,'FaceColor','r','FaceAlpha',.2);
    
    save('..\data\prostate\prostate_034_w_obj.mat', 'faces', 'vertices' );
    Objs = {'..\data\prostate\prostate_034_w_obj.mat'};
    
    outName_parameterization = SpharmMatParameterization(paramCALD, Objs, 'ParamCALD');
    outNames_expansion = SpharmMatExpansion(paramCALD, outName_parameterization, 'ExpLSF');
    load( '..\data\prostate\prostate_034_w_CALD_LSF_des.mat' );
    
    
    metric_w = [];
    for alpha=-pi/2:step_rot:pi/2
        rvec = spharm_rotate(alpha, 0.0, 0.0, full_fvec, paramCALD.MaxSPHARMDegree);
        metric_w = [metric_w; norm(fvec - rvec)];
    end
    alpha=-pi/2:step_rot:pi/2;
    h = figure;
    
    set(gca,'YGrid','on',...
        'XTickLabel',{'-pi/2','-pi/4','0','pi/4','pi/2'},...
        'XTick',[-1.571 -0.785 0 0.785 1.571],...
        'XGrid','on',...
        'FontSize',20);
    
    % xlim(axes1,[-1.571 1.571]);
    
    % Create plot
    plot(alpha,metric_w,'LineWidth',2), title(['w = ' num2str(w)]);
    
    % Create title
    %title('No Missing Data','FontSize',20);
    
    % Create xlabel
    xlabel('Rotation about the lateral axis of the prostate in radians',...
        'FontSize',20);
    
    % Create ylabel
    ylabel('SPHARM distance metric','FontSize',20);
    close(h);
end
