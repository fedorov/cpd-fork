add_bcpd_paths;
clear; close all; clc;

%% Load the sphere
[p, t] = read_obj('..\data\synthetic\sphere.obj');
d=0.4;

%% Calculate the radius of the sphere and its center
c = mean(p);
pp = bsxfun(@minus, p, c);
r = sqrt(trace((pp')*pp)/size(p,1));

%% Calculate the volume of intersection
% See http://mathworld.wolfram.com/Sphere-SphereIntersection.html
V_int = pi*(2*r-d)^2;
V_int = V_int*(d^2+4*d*r);
V_int = V_int/(12*d);

V_sphere = 4*pi*r^3/3;

dice_ground = 2*V_int/(2*V_sphere);

%% Visualize the sphere and translate it by e_s units
h(1) = patch('Vertices',p,'Faces',t,'FaceColor','r','FaceAlpha',.2);
pp = bsxfun(@plus, p, [d 0 0]);
h(2) = patch('Vertices',pp, 'Faces',t,'FaceColor','b','FaceAlpha',.2);
axis equal vis3d; drawnow;

%% Calculate dice using the "inpolyhedron" technique
tic;
dice_polyhedron = cpd_dice(p,t,pp,t,50,50,50);
dice_time = toc;
fprintf('InPolyhedron dice estimate: time=%f, val=%f\n',dice_time, dice_polyhedron); 

%% Calculate dice using CSG
tic;
[dice_csg, vint, v1, v2] = dice(p', t', pp', t', 1e-7);
dice_time = toc;
fprintf('CSG dice estimate: time=%f, val=%f\n',dice_time, dice_csg);

%% Calculate dice using VolumeOp binary
tic;
V_1 = abs(VolumeSurface(p,t));
V_2 = abs(VolumeSurface(pp,t));
% Write out the two patches into off files
write_off('..\data\sphere1.off',p,t); 
write_off('..\data\sphere2.off',pp,t);
% Calculate the intersection
system('..\Core\VolumeOp.exe ..\data\sphere1.off ..\data\sphere2.off ..\data\int_sphere.off');
% Read in the intersection
[pint,fint] = read_off('..\data\int_sphere.off');
write_obj(pint', fint', '..\data\int_sphere.obj');
V_int = abs(VolumeSurface(pint',fint'));
% Calculate Dice
dice_volumeop = 2.0*V_int/(V_1+V_2);
dice_time = toc;
fprintf('VolumeOp dice estimate: time=%f, val=%f\n',dice_time, dice_volumeop);