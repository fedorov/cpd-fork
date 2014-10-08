
[Y,fY] = read_ply('P30\US_affine_model.ply');
[X,fX] = read_ply('P30\MR_model.ply');

numberOfFaces = 2000;

pX.faces = fX;
pX.vertices = X;

pY.faces = fY;
pY.vertices = Y;

nfx = reducepatch(pX,numberOfFaces);
nfy = reducepatch(pY,numberOfFaces);

% [TY, ~, ~, ~, ~, newSigma2, ~, fem, u, ~] = cpd_fem_only(nfx.vertices, nfy.vertices, nfy.faces, w, errtol, maxiters, eye(3), [0;0;0], 1.0, sigma2, beta, E, nu, [], [], []); 

X = nfx.vertices;
XX = cut_data(X, 0.2, 0.2, [0 0 1]);
clf;
plot3(X(:,1), X(:,2), X(:,3),'ob','MarkerSize', 5, 'MarkerFaceColor', 'b');
hold on;
plot3(XX(:,1), XX(:,2), XX(:,3),'or','MarkerSize', 5, 'MarkerFaceColor', 'r');
hold off;