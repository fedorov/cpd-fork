%% Tests registration of a set of points using a simple SVD technique
add_bcpd_paths

%% Generate points

N = 100;

R = rotation_matrix(pi/6,pi/6,0);
X = rand(N, 3)*diag([10 5 3])*(R)+1;

R2 = eye(3); %rotation_matrix(rand(1)*2*pi, rand(1)*2*pi,rand(1)*2*pi);
t = [1 2 3]; % rand(1,3);
s = 1.5; % 1-rand(1)*0.5;
Y = bsxfun(@plus, s*X*R2, t);

%% Register
[R, t, s, TY] = tight_box_reg(X, Y);


%% Plot
plot3(X(:,1), X(:,2), X(:,3), 'b.', 'MarkerSize', 5);
hold on;
plot3(Y(:,1), Y(:,2), Y(:,3), 'k.', 'MarkerSize', 5);
plot3(TY(:,1), TY(:,2), TY(:,3), 'r.', 'MarkerSize', 20);
legend('X', 'Y', 'TY');

% show boxes
%% box
[R1, c1, w1, B1] = tight_box(X);
[R2, c2, w2, B2] = tight_box(TY);

% patch faces
F = [ 1 2 3 4;
      5 8 7 6;
      1 5 6 2;
      4 3 7 8;
      1 4 8 5;
      2 6 7 3];
box.faces = F;
box.vertices = B1;
patch(box, 'FaceColor', 'blue', 'FaceAlpha', 0.1);
box.vertices = B2;
patch(box, 'FaceColor', 'red', 'FaceAlpha', 0.1);
axis equal;
hold off;



