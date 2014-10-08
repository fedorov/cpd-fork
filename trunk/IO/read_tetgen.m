function [ nodes, elems ] = read_tetgen( nodefile, elefile )
%READ_TETGEN Reads data from tetgen .ele and .node files
%
%   [nodes, elems] = read_tetgen(nodefile, elefile)
%       Parses tetgen output files for creating a tetrahedral
%       finite element model.
%
%       nodefile:   tetgen .node file
%       elefile:    tetgen .ele file
%       nodes:      Nx3 array of node locations
%       elems:      Mx4 array of node indices that make up elements
%
%   Copyright 2013 C. Antonio Sanchez [antonios@ece.ubc.ca]

%% Load node file
% First line: 
%   <# points> <dimension=3> <# attributes> <# boundary markers (0 or 1)>
% Remaining lines list # of points:
%   <point #> <x> <y> <z> [attributes] [boundary marker]

fnode = fopen(nodefile,'r');
if (fnode < 0)
    error(['Cannot open the file ''', nodefile,'''']);
end

str = fgetl(fnode);
A = sscanf(str,'%d %d %d %d');

N = A(1);
idxs = zeros(N,1);  % use idxs to resolve nodes
nodes = zeros(N,3);
for i=1:N
    str = fgetl(fnode);
    A = sscanf(str, '%d %f %f %f');
    idxs(i) = A(1);
    nodes(i,:) = A(2:4);
end
fclose(fnode);

%% Load elem file
% First line: 
%   <# tetrahedra> <nodes per tet> <# of attributes>
% Remaining lines list of # of tetrahedra:
%   <tetrahedron #> <node> <node> <node> <node> ... [attributes]

felem = fopen(elefile, 'r');
if (felem < 0)
    error(['Cannot open the file ''', elefile,'''']);
end

str = fgetl(felem);
A = sscanf(str, '%d %d %d');
M = A(1);
D = A(2);

fmtStr = ['%d', repmat(' %f',[1, D]) ];

elems = zeros(M,D);
for i=1:M
    str = fgetl(felem);
    A = sscanf(str, fmtStr);
    A = A(2:(D+1));
    
    % replace indices based on index map
    [~, e] = ismember(A,idxs);
    elems(i,:) = e(:);
end

fclose(felem);

end

