% VOLUMESURFACE            Finds the volume inside a triangular surface mesh
%
% V = VolumeSurface(VER,FACE)
%
% Finds the volume of a watertight surface formed by the vertices VER (Px3)
% connected by the faces FACE (Tx3). All face normals must direct the same
% side of the volume (you get a +/- result depending on this)
%
% The divergence theorem is used, to infer the volume from the surface
% integral of a uniformly-increasing vector function [e.g., F(x,y,z)={x,0,0}].
% The triangles are discretized at their centroids, to compute the "flow"
% of this function across the surface (their x cross-sections, indeed).
%
% Copyright (c) Orcun Goksel - Nov 2008
%

function vol = VolumeSurface(p,t)

% Form a 3D array of triangle corners
v = reshape(p(t',:)',3,3,[]);

% Find face normals, with magnitude twice the areas
nc = cross(diff(v(:,[1 2],:),1,2),diff(v(:,:,:),2,2));

% cross-section along x is multiplied by barycenter function value
vol = squeeze(nc(1,1,:))' * squeeze(mean(v(1,:,:),2)) /2;

end