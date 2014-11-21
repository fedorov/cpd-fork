%% Tight-fitting box test

add_bcpd_paths

%% Generate points

N = 100;
R = rotation_matrix(pi/6,pi/6,0);
X = rand(N, 3)*diag([10 5 3])*(R')+1;

%% box
[R, c, w, B] = tight_box(X);

% patch faces
F = [ 1 2 3 4;
      5 8 7 6;
      1 5 6 2;
      4 3 7 8;
      1 4 8 5;
      2 6 7 3];
      
labels = {'1','2','3','4','5','6','7','8'};
      
clf;
plot3(X(:,1), X(:,2), X(:,3),'.b','MarkerSize',10);
hold on;
plot3(B(:,1), B(:,2), B(:,3),'.r','MarkerSize',20);
text(B(:,1)+0.1, B(:,2)+0.1, B(:,3)+0.1, labels);

box.vertices = B;
box.faces = F;
patch(box, 'FaceColor', 'red', 'FaceAlpha', 0.1);
axis equal;
hold off;