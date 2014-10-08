% Written by Abtin Rasoulian: abtinr@ece.ubc.ca
function writePolyDataInVTK(vertices, faces, fileName)
fid = fopen(fileName, 'w');
fprintf(fid, '# vtk DataFile Version 3.0\nvtk output \nASCII \nDATASET POLYDATA \nPOINTS %d float\n', size(vertices, 1));
fprintf(fid, '%2.4f %2.4f %2.4f\n', vertices');
fprintf(fid, '\nPOLYGONS %d %d \n', size(faces, 1), size(faces, 1)*4);
D = [ones(size(faces, 1), 1)*3 faces-1];
fprintf(fid, '%d %d %d %d\n', D');
fprintf(fid, '\nCELL_DATA %d \nPOINT_DATA %d', size(faces, 1), size(vertices, 1));
fclose(fid);

