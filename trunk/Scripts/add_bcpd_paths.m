p_generated = genpath([root '/Core']);
addpath(p_generated);
addpath([root '/IO']);
addpath([root '/data']);
addpath([root '/Scripts']);

p_generated = genpath([root '/ThirdParty/tetgen1.4.3']);
addpath(p_generated);

p_generated = genpath([root '/ThirdParty/maslib']);
addpath(p_generated);

p_generated = genpath([root '/ThirdParty/spharm']);
addpath(p_generated);

addpath([root '/ThirdParty/kdtree']);
addpath([root '/ThirdParty/rbf']);

clear p_generated;
