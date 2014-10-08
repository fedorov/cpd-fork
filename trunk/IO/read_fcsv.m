function P = read_fcsv( filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(filename,'r','b');

fgetl(fid); % # Markups fiducial file version = 4.3
fgetl(fid); % CoordinateSystem = 0
fgetl(fid); % columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID

% Read the first line
% str = fgetl(fid);
% if (str(1,30) == ',')
%     first_fiducial = str2num(str(1,28:29));
% elseif (str(1,29) == ',')
%     first_fiducial = str2num(str(1,28));
% elseif (str(1,31) == ',')
%     first_fiducial = str2num(str(1,28:30));
% else
%     error('Something unexpected happened. Please check the fiducial file.');
% end

% string = ['vtkMRMLMarkupsFiducialNode_' num2str(first_fiducial) ',%f,%f,%f'];
% P(i,:) = sscanf(str, string);
P = [];
i = 1; 
while ~feof(fid)
    string = fgetl(fid);
    comma = strfind(string, ',');
    string(1:comma) = [];
    P(i,:) = sscanf(string, '%f,%f,%f');  %
    i=i+1;
end

fclose(fid);
end

