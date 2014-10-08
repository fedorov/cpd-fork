% Author: Abtin Rasoulian: abtinr@ece.ubc.ca
function [origin spacing dimension] = getVTKFileInfo(fileName)
fid = fopen(fileName, 'rt');
for i = 1:4
    fgetl(fid);
end
for i=1:3
    str = fgetl(fid);
    if(findstr('DIMENSION', str))
        dimension = str2num(str(11:length(str)));
    elseif(findstr('ORIGIN', str))
        origin = str2num(str(8:length(str)));
    elseif(findstr('SPACING', str))
        spacing = str2num(str(8:length(str)));
    end
end

fclose(fid);
