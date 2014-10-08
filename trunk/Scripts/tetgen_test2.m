add_bcpd_paths;

%% Load the SSM and write out the mean shape to a .ply file
load('..\data\SSM.mat');
write_ply(SSM.mean, SSM.faces, '..\data\mean.ply');

%% Call tetgen to create a tetrahedral mesh
system( '..\ThirdParty\tetgen1.4.3\build\Release\tetgen.exe -q2 -p ..\data\mean.ply -O');

%% Read the constructed mesh
[v,f] = read_off('..\data\mean.1.off');
patch('Vertices',v','Faces',f','FaceColor','r','FaceAlpha',.5);