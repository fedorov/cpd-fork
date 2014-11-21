add_bcpd_paths;
% matlabpool open;

%% load data
[X.vertices,X.faces] = readPolyDataInVTK('..\data\prostate\P1.vtk');

%% create model
fem = createBoundingFem(X.vertices, 0.5, [5, 5, 5]);
nodes = getNodes(fem);

figure(1);
hold off;
plot3(nodes(:,1), nodes(:,2), nodes(:,3),'.r','MarkerSize',10);
hold on;
patch(X,'FaceColor','blue','FaceAlpha',0.2, 'EdgeColor', 'blue', 'EdgeAlpha',0.2);
axis equal;
drawnow;

EIn = 10;
EOut = 5;
nuIn = 0.49;
nuOut = 0.49;

mat = fem_material_spatial_linear_twomaterials(EIn, nuIn, EOut, nuOut, X.vertices, X.faces);
setMaterial(fem, mat);

%% Containing element


Phi = getInterpolationMatrix(fem, X.vertices);
% K = getStiffnessMatrix(fem);




%% cleanup
delete(fem);