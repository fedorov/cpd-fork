add_bcpd_paths;

%% Run tetgen.exe with example.poly to ensure correct installation
system( '..\ThirdParty\tetgen1.4.3\build\Release\tetgen.exe -q2 -C ..\data\example.poly' );

%% This should produce the following files in your data folder:
% example.1.ele example.1.face example.1.node