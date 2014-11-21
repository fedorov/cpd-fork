add_bcpd_paths;

%% Load data
load('..\data\SSM.mat');
[X,Xf] = readPolyDataInVTK('..\data\prostate\P1.vtk');
Y = SSM.mean;
Yf = SSM.faces;

FV.faces = Xf;
FV.vertices = X;

[az, el] = view;
clf;
patch(FV,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
hold on;
axis([-3,3,-3,3,-3,3]);
% without deformation
plot3(Y(:,1), Y(:,2), Y(:,3),'ob','MarkerSize', 5, 'MarkerFaceColor', 'b');
axis([-3,3,-3,3,-3,3]);
legend({'target','ssm'},'location','NEo')
view(az, el);
drawnow;

%% CPD Parameters
w = 0.01;
errtol = 1e-10;
maxiters = 100;

lambda = 0.1;
beta2 = 3.5;
sigma2 = [];


%% FEM Parameters
E = 480;
nu = 0.3;

%% CPD + Fem from boundary
[TY, P, sigma2, fem, U, Phi] = cpd_plus_boundary_fem(X, Y, Yf, lambda, beta2, w, errtol, maxiters, sigma2, E, nu );

%% plot CPD results
clf;
subplot(1,2,1);
patch(FV,'FaceColor','red','FaceAlpha',0.2, 'EdgeColor', 'red', 'EdgeAlpha',0.2);
hold on;
axis([-3,3,-3,3,-3,3]);
% without deformation
plot3(TY(:,1), TY(:,2), TY(:,3),'ob','MarkerSize', 5, 'MarkerFaceColor', 'b');
axis([-3,3,-3,3,-3,3]);
legend({'target','ssm'},'location','NEo')
view(az, el);
drawnow;

%% Interpolate some points

% original FEM
setSurfaceMesh(fem, Yf);

% new FEM with moved node locations (for reverse interpolation)
nodes2 = getNodes(fem);
elems2 = getElements(fem);
nodes2 = nodes2 + U;
femNew = fem_model(nodes2, elems2);
setSurfaceMesh(femNew, Yf);

% pick some random points from first model
Npoints = 10;
[ R, c, w, B ] = tight_box( Y );
Ypnts = rand(Npoints, 3);
Ypnts = (Ypnts - 0.5).*repmat(w, [size(Ypnts,1) 1]);
Ypnts = Ypnts*R + repmat(c, [size(Ypnts,1) 1]);
[P, in]= getInterpolationMatrix(fem, Ypnts);

% find locations in second model
Xpnts = P*getNodes(femNew);

% Verify correctness
[P2, in2] = getInterpolationMatrix(femNew, Xpnts);
Ypnts2 = P2*getNodes(fem);

% plot:
subplot(1,2,2);
plot3(Ypnts(:,1), Ypnts(:,2), Ypnts(:,3),'ob','MarkerSize', 5, 'MarkerFaceColor', 'b');
hold on;
axis([-3,3,-3,3,-3,3]);
% interpolated
plot3(Xpnts(:,1), Xpnts(:,2), Xpnts(:,3),'or','MarkerSize', 5, 'MarkerFaceColor', 'r');
axis([-3,3,-3,3,-3,3]);

legend({'Y Points','X points'},'location','NEo')
view(az, el);
drawnow;

% Conclusion: works for interior points