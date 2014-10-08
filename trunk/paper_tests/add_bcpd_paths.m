p_generated = genpath('../Core');
addpath(p_generated);
addpath('../IO');
addpath('../data');
addpath('../Scripts');

p_generated = genpath('../ThirdParty/tetgen1.4.3');
addpath(p_generated);

p_generated = genpath('../ThirdParty/maslib');
addpath(p_generated);

addpath('../ThirdParty/kdtree');
addpath('../ThirdParty/rbf');

clear p_generated;