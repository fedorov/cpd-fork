% Author: Abtin Rasoulian, abtinr@ece.ubc.ca
function m_fv = extractSurface(fileName, isovalue)
if(nargin == 1)
    isovalue = 0.5;
end
% isovalue = 1200;
V = readVTK(fileName);
% V(V > 64000) = 0;%65535) = 0;
[origin spacing dimension] = getVTKFileInfo(fileName);
[X Y Z] = meshgrid( ...
    origin(2):spacing(2):origin(2)+spacing(2)*(size(V,2)-1), ...
    origin(1):spacing(1):origin(1)+spacing(1)*(size(V,1)-1), ...
    origin(3):spacing(3):origin(3)+spacing(3)*(size(V,3)-1));

fv = isosurface(X,Y,Z,V,isovalue,'verbose');

m_fv = fv;
m_fv.vertices(:,1) = fv.vertices(:,2);
m_fv.vertices(:,2) = fv.vertices(:,1);
% m_fv = m_fv.vertices;