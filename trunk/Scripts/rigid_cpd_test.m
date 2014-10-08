add_bcpd_paths;

%% Test values
N = 50;
M = 45;
D = 3;

% random points
X = rand(N,D);

% random rotation
alpha = 1;
[R, T] = qr(randn(D,D)+alpha*eye(D));
R = R*diag(sign(diag(T)));
if (det(R) < 0)
    R(:,1) = -R(:,1);
end
% R = eye(D);

% random scale [-1.5 1.5]
% s = 1;
s = 2*((2*rand(1)-1)*0.5+1);


% random translation
% t = zeros(D,1);
t = rand(D,1);


% choose random set of X, then transform
Yidx = randperm(N,M);
Y = X(Yidx,:);
Y = (Y - repmat(t', [M,1]))*R/s;

%% CPD values
w = 0.01;
errtol = 1e-15;
maxiters = 100;
sigma2 = 0.01; %[]; %0.01;

%% Test rigid cpd
[TY, Rout, tout, sout, Pout] = cpd_rigid(X, Y, w, errtol, maxiters, [],[],[],sigma2);
[TY2, Rout2, tout2, sout2, Pout2] = cpd_rigid_symmetric(X, Y, w, errtol, maxiters, [],[],[],sigma2);

subplot(1,2,1);
plot3(X(:,1), X(:,2), X(:,3), 'b.', 'MarkerSize', 5);
hold on;
plot3(Y(:,1), Y(:,2), Y(:,3), 'k.', 'MarkerSize', 5);
plot3(TY(:,1), TY(:,2), TY(:,3), 'r.', 'MarkerSize', 20);
legend('X', 'Y', 'TY');
title('Rigid CPD');
hold off;

subplot(1,2,2);
plot3(X(:,1), X(:,2), X(:,3), 'b.', 'MarkerSize', 5);
hold on;
plot3(Y(:,1), Y(:,2), Y(:,3), 'k.', 'MarkerSize', 5);
plot3(TY2(:,1), TY2(:,2), TY2(:,3), 'r.', 'MarkerSize', 20);
legend('X', 'Y', 'TY');
title('Symmetric rigid CPD');
hold off;

disp('[ Rotation | Estimated Rotation ]');
disp( [R Rout] );

disp('[ Translation | Estimated Translation ]');
disp( [t tout] );

disp('[ Scale | Estimated Scale ]');
disp( [s sout] );

disp('[ Alignment Index | Estimated Alignment ]');
[~, Yidxout] = max(Pout,[],2);

disp( [Yidx' Yidxout] );