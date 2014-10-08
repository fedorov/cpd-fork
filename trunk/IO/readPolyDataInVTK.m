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
i=1;
while( le(i,pointNr) )
    str=fgetl(fid);
    temp=str2num(str);
    if ( length(temp) == 3 )
        vertices(i,:) = temp;
        i=i+1;
    elseif( length(temp) == 6 )
        vertices(i,:) = [temp(1) temp(2) temp(3)];
        i=i+1;
        vertices(i,:) = [temp(4) temp(5) temp(6)];
        i=i+1;
    elseif( length(temp) == 9 )
        vertices(i,:) = [temp(1) temp(2) temp(3)];
        i=i+1;
        vertices(i,:) = [temp(4) temp(5) temp(6)];
        i=i+1;
        vertices(i,:) = [temp(7) temp(8) temp(9)];
        i=i+1;
    else
         error('Something went wrong while reading vertices');
    end
end

s = fgetl(fid);
if (isempty(s))
    s = fgetl(fid);
end
faceNr = sscanf(s, '%*s%d');
faces = zeros(faceNr, 3);
for i=1:faceNr
    tmp = str2num(fgetl(fid));
    faces(i,:) = tmp(2:4)+1;
end
fclose(fid);

