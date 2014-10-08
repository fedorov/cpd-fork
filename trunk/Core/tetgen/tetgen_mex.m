% TETGEN Interface to the tetgen library.
% 
% [nodes,faces,tets] = tetgen(V, F, quality, switches)
%   V:  3xN set of vertices
%   F:  KxM set of faces, where K is the number of vertices per face
%       (optional, defaults to computing convex hull)
%   quality:  quality of tet mesh (optional, default = 2)
%   switches: tetgen switches     (optional, default = 'Q')
%   nodes: 3xNn set of nodes
%   faces: 3xNf set of triangular faces (or convex hull if F is empty)
%   tets:  4xMt set of tetrahedral elements
