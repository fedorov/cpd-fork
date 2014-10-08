add_bcpd_paths;

%% Test values
N = 10;
M = 10;
D = 3;

% random points
X = rand(N,D);

% random affine
alpha = 1;
A = randn(D,D)+alpha*eye(D);
% A = eye(D);

% random translation
% t = zeros(D,1);
t = rand(D,1);


% choose random set of X, then transform
Yidx = randperm(N,M);
Y = X(Yidx,:);
Y = (Y - repmat(t', [M,1]))/(A');

%% CPD values
w = 0.0;
errtol = 1e-15;
maxiters = 100;
sigma2 = []; %0.01;

%% Test affine cpd
[TY, Aout, tout, Pout] = cpd_affine(X, Y, w, errtol, maxiters, [],[],sigma2);

plot3(X(:,1), X(:,2), X(:,3), 'b.', 'MarkerSize', 5);
hold on;
plot3(Y(:,1), Y(:,2), Y(:,3), 'k.', 'MarkerSize', 5);
plot3(TY(:,1), TY(:,2), TY(:,3), 'r.', 'MarkerSize', 20);
legend('X', 'Y', 'TY');
hold off;

disp('[ Affine | Estimated Affine ]');
disp( [A Aout] );

disp('[ Translation | Estimated Translation ]');
disp( [t tout] );

disp('[ Alignment Index | Estimated Alignment ]');
[~, Yidxout] = max(Pout,[],2);

disp( [Yidx' Yidxout] );