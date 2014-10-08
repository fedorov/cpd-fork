add_bcpd_paths;

%% Prefix and suffixes
prefix = 'P';
segSuffix = 'Seg.vtk';
volSuffix = 'Volume.vtk';
mr = 'MR';
us = 'US';

%% Registration parameters
s=1; R=eye(3); t=[0;0;0];
w=0.0; errtol=1e-10; maxiters=100; sigma2=[];
beta=1.0; E=4.80; nu=0.49;

%% Data folder and file names
directory = '..\data\prostate\';
flist = fopen('..\data\prostate\MR_TRUS.txt','r');

str='';
while( ~strcmp(str,'Folders') )
    str=fgetl(flist);
end

noPatients = str2num(fgetl(flist));
dPatients = zeros(noPatients,3);

for i=1:noPatients
    dPatients(i,:) = fgetl(flist);
end

fid = fopen('..\data\prostate\MR_TRUS.log','wt');
formatSpec = 'Processing Patient %s\n';

% First column before any registration
% Second column rigid only and third rigid + Deformable
DiceArray = zeros(noPatients,3); 

%% Read the data and perform registration
for i=1:noPatients
    fprintf(fid,formatSpec,dPatients(i,:));
    
    [X,fX] = readPolyDataInVTK([directory prefix dPatients(i,:) '_' us '_' segSuffix]);
    [Y,fY] = readPolyDataInVTK([directory prefix dPatients(i,:) '_' mr '_' segSuffix]);
    
    % Read the MR image to resample in TRUS
    [MR,origin,spacing] = readImageDataInVTK([directory prefix dPatients(i,:) '_' mr '_' volSuffix]);
    
    % Calculate dice before registration
    VX = abs(VolumeSurface(X,fX));
    VY = abs(VolumeSurface(Y,fY));
    write_off('..\data\prostate\X.off',X,fX);
    write_off('..\data\prostate\Y.off',Y,fY);
    
    system('..\Core\VolumeOp.exe ..\data\prostate\X.off ..\data\prostate\Y.off ..\data\prostate\int_XY.off');
    [pint,fint] = read_off('..\data\prostate\int_XY.off');
    VXY = abs(VolumeSurface(pint',fint'));
    
    dice_before = 2*VXY/(VX+VY);
    
    % Perform rigid only registration
    tic;
    [TY, ~, ~, ~, ~] = cpd_rigid(X, Y, w, errtol, maxiters, [],[],0,sigma2);
    time=toc; fprintf('Rigid registration time is %f seconds\n', time);
    
    % Calculate dice after registration
    write_off('..\data\prostate\TY.off',TY,fY);
    system('..\Core\VolumeOp.exe ..\data\prostate\X.off ..\data\prostate\TY.off ..\data\prostate\int_XTY.off');
    [pint,fint] = read_off('..\data\prostate\int_XTY.off');
    VTY = abs(VolumeSurface(TY,fY));
    VXTY = abs(VolumeSurface(pint',fint'));
    
    dice_rigid = 2*VXTY/(VX+VTY);
    
    % Perform rigid + deformable registration
    tic;
    [TY, ~, ~, ~, fem, u, outsigma, ~] = cpd_rigid_fem(X, fX, Y, fY, w, errtol, maxiters, R, t, s, sigma2, beta, E, nu);
    time=toc; fprintf('FEM registration time is %f seconds\n', time);
    
%     maxY = max(Y,[],1);
%     minY = min(Y,[],1);
%     width = (maxY-minY);
%     
%     rx=10; ry=10; rz=10;
%     
%     sp = [2*width(1)/(rx-1), 2*width(2)/(ry-1), 2*width(3)/(rz-1)];
%     x1 = -width(1):sp(1):width(1);
%     x2 = -width(2):sp(2):width(2);
%     x3 = -width(3):sp(3):width(3);
%     
%     [x1g,x2g,x3g] = meshgrid(x1, x2, x3);
%     xg = [reshape(x1g,1,[])' reshape(x2g,1,[])' reshape(x3g,1,[])'];
%     
%     tic;
%     phi = getInterpolationMatrix(fem, xg);
%     time=toc; fprintf('Interpolation time is %f seconds\n', time);
%     
%     dg = phi*u;
    
    % Calculate dice after registration
    write_off('..\data\prostate\TY.off',TY,fY);
    system('..\Core\VolumeOp.exe ..\data\prostate\X.off ..\data\prostate\TY.off ..\data\prostate\int_XTY.off');
    [pint,fint] = read_off('..\data\prostate\int_XTY.off');
    VTY = abs(VolumeSurface(TY,fY));
    VXTY = abs(VolumeSurface(pint',fint'));
    
    dice_fem = 2*VXTY/(VX+VTY);
    
    DiceArray(i,1) = dice_before;
    DiceArray(i,2) = dice_rigid;
    DiceArray(i,3) = dice_fem;
end

save('..\data\prostate\DiceFEM.mat','DiceArray');
fclose(fid);