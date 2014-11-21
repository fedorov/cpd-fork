add_bcpd_paths;

%% FISH

fish = [-0.91542, -0.16535, 0.00000;
       -0.89051, -0.10213, 0.00000;
       -0.85730, -0.10213, 0.00000;
       -0.82408, -0.14955, 0.00000;
       -0.87390, -0.19696, 0.00000;
       -0.92372, -0.16535, 0.00000;
       -1.06488, -0.29180, 0.00000;
       -1.02336, -0.16535, 0.00000;
       -0.98185, -0.05471, 0.00000;
       -0.95694, 0.02432, 0.00000;
       -0.91542, 0.11915, 0.00000;
       -0.88221, 0.18237, 0.00000;
       -0.80747, 0.26140, 0.00000;
       -0.75765, 0.34043, 0.00000;
       -0.67462, 0.46687, 0.00000;
       -0.58328, 0.54590, 0.00000;
       -0.51686, 0.60912, 0.00000;
       -0.40061, 0.70396, 0.00000;
       -0.33418, 0.81460, 0.00000;
       -0.26775, 1.16232, 0.00000;
       -0.20133, 1.46262, 0.00000;
       -0.15151, 1.76293, 0.00000;
       -0.14320, 1.88938, 0.00000;
       -0.11829, 2.11065, 0.00000;
       -0.10169, 1.96840, 0.00000;
       -0.09338, 1.74713, 0.00000;
       -0.09338, 1.54165, 0.00000;
       -0.08508, 1.30457, 0.00000;
       -0.09338, 1.11490, 0.00000;
       -0.07678, 0.90943, 0.00000;
       -0.07678, 0.71976, 0.00000;
       -0.04356, 0.56171, 0.00000;
       -0.03526, 0.48268, 0.00000;
       0.04777, 0.40365, 0.00000;
       0.11420, 0.34043, 0.00000;
       0.16402, 0.27721, 0.00000;
       0.22215, 0.21398, 0.00000;
       0.29688, 0.10334, 0.00000;
       0.36330, -0.00729, 0.00000;
       0.44634, -0.08632, 0.00000;
       0.50446, -0.05471, 0.00000;
       0.57089, 0.04012, 0.00000;
       0.60410, 0.18237, 0.00000;
       0.67053, 0.38784, 0.00000;
       0.72865, 0.62493, 0.00000;
       0.81999, 0.67235, 0.00000;
       0.90302, 0.76718, 0.00000;
       0.88642, 0.62493, 0.00000;
       0.86981, 0.46687, 0.00000;
       0.86151, 0.26140, 0.00000;
       0.83660, 0.10334, 0.00000;
       0.84490, -0.13374, 0.00000;
       0.84490, -0.30760, 0.00000;
       0.84490, -0.52888, 0.00000;
       0.87811, -0.73435, 0.00000;
       0.92794, -0.89241, 0.00000;
       0.98606, -1.00305, 0.00000;
       1.03588, -1.12949, 0.00000;
       1.04418, -1.20852, 0.00000;
       0.96945, -1.14530, 0.00000;
       0.88642, -1.05046, 0.00000;
       0.79508, -0.97144, 0.00000;
       0.68714, -0.92402, 0.00000;
       0.59580, -0.82919, 0.00000;
       0.53768, -0.75016, 0.00000;
       0.47125, -0.63952, 0.00000;
       0.41312, -0.57630, 0.00000;
       0.28857, -0.52888, 0.00000;
       0.18893, -0.46566, 0.00000;
       0.08929, -0.46566, 0.00000;
       -0.03526, -0.49727, 0.00000;
       -0.08508, -0.49727, 0.00000;
       -0.00205, -0.57630, 0.00000;
       0.09760, -0.65532, 0.00000;
       0.17233, -0.68694, 0.00000;
       0.26366, -0.73435, 0.00000;
       0.21384, -0.73435, 0.00000;
       0.03947, -0.73435, 0.00000;
       -0.15151, -0.73435, 0.00000;
       -0.26775, -0.71855, 0.00000;
       -0.40061, -0.67113, 0.00000;
       -0.46704, -0.60791, 0.00000;
       -0.55837, -0.52888, 0.00000;
       -0.59989, -0.51307, 0.00000;
       -0.77426, -0.46566, 0.00000;
       -0.89051, -0.48146, 0.00000;
       -0.98185, -0.43405, 0.00000;
       -1.03997, -0.37082, 0.00000;
       -1.05658, -0.29180, 0.00000;
       -0.05048, -0.75806, 0.00000;
       0.09967, -0.75806, 0.00000];
%% Bounding box for fish and control points
up_bound = max(fish);
low_bound = min(fish);

numberNodes = 10;

[ctX,ctY,ctZ] = meshgrid(low_bound(1):(up_bound(1)-low_bound(1))/numberNodes:up_bound(1), low_bound(2):(up_bound(2)-low_bound(2))/numberNodes:up_bound(2), low_bound(3):(up_bound(3)-low_bound(3))/numberNodes:up_bound(3));
   
%% Transform parameters
s = 1.0;
R = eye(3);
t = [0 0 1]';

%% Inputs
N = size(fish,1);
M = N;

Yidx = (1:M)';

X = s*fish*(R')+repmat(t', [N, 1]);
Y = fish(Yidx,:);

% Registration parameters
w = 0.05;
errtol = 1e-3;
maxiters = 100;
lambda = 0.001;
sigma2 = [];

%% Run cpd
[ TY, P, sigma2 ] = tps_rpm( X, Y, lambda, w, errtol, maxiters, sigma2)

%% Display results
figure, plot3(X(:,1), X(:,2), X(:,3), 'b.', 'MarkerSize', 5);
hold on;
plot3(Y(:,1), Y(:,2), Y(:,3), 'g.', 'MarkerSize', 5);
plot3(TY(:,1), TY(:,2), TY(:,3), 'r.', 'MarkerSize', 20);
legend('X', 'Y', 'TY', 'location', 'SO');
title('TPS-RPM');
hold off;