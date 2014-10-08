% Written by Abtin Rasoulian: abtinr@ece.ubc.ca
function [vertices, faces] = readPolyDataInVTK(fileName)
fid = fopen(fileName, 'r','b');
fgetl(fid); % # vtk DataFile Version x.x
fgetl(fid);
fgetl(fid);
fgetl(fid);

s = fgetl(fid);
pointNr = sscanf(s, '%*s%d');
vertices = zeros(pointNr, 3);
for i=1:pointNr
    vertices(i,:) = str2num(fgetl(fid));
end
fgetl(fid);
s = fgetl(fid);
faceNr = sscanf(s, '%*s%d');
faces = zeros(faceNr, 3);
for i=1:faceNr
    tmp = str2num(fgetl(fid));
    faces(i,:) = tmp(2:4)+1;
end

