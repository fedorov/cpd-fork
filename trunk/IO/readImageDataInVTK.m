function [V,origin,spacing] = readImageDataInVTK(vtkfile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: V = readVTK(vtkfile)
%
%   V:       The matrix to be stored
%   vtkfile: The filename
%   notes:   Only reads binary STRUCTURED_POINTS
%
% Erik Vidholm 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = 0;

% open file (OBS! big endian format)
fid = fopen(vtkfile,'r','b');

if( fid == -1 )
  return
end

fgetl(fid); % # vtk DataFile Version x.x
fgetl(fid); % comments
fgetl(fid); % BINARY
fgetl(fid); % DATASET STRUCTURED_POINTS

s = fgetl(fid); % DIMENSIONS NX NY NZ
sz = sscanf(s, '%*s%d%d%d').';

str = fgetl(fid);
spacing = sscanf(str, 'SPACING %f %f %f'); % SPACING SX SY SZ

if (isempty(spacing))
origin = sscanf(str, 'ORIGIN %f %f %f');  % ORIGIN OX OY OZ
spacing = sscanf(fgetl(fid), 'SPACING %f %f %f'); % SPACING SX SY SZ
else
    origin = sscanf(fgetl(fid), 'ORIGIN %f %f %f');  % ORIGIN OX OY OZ
end

fgetl(fid); % POINT_DATA NXNYNZ

s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
svstr = sscanf(s, '%s', 1);
dtstr = sscanf(s, '%*s%*s%s');
dtstr_vec = sscanf(s, '%*s%*s%*s%s');

if( strcmp(svstr,'SCALARS') > 0 )
  fgetl(fid); % the lookup table
  if( strcmp(dtstr,'unsigned_char') > 0 ) 
    % read data
    V = fread(fid,prod(sz),'*uint8');
    V = reshape(V,sz);
  elseif( strcmp(dtstr,'char') > 0 )
    % read data
    V = fread(fid,prod(sz),'*int8');
    V = reshape(V,sz);
  elseif( strcmp(dtstr,'unsigned_short') > 0 )
    % read data
    V = fread(fid,prod(sz),'*uint16');
    V = reshape(V,sz);
  elseif( strcmp(dtstr,'short') > 0 )
    % read data
    V = fread(fid,prod(sz),'*int16');
    V = reshape(V,sz);
  elseif( strcmp(dtstr,'unsigned_int') > 0 )
    % read data
    V = fread(fid,prod(sz),'*uint32');
    V = reshape(V,sz);
  elseif( strcmp(dtstr,'int') > 0 )
    % read data
    V = fread(fid,prod(sz),'*int32');
    V = reshape(V,sz);
  elseif( strcmp(dtstr,'float') > 0 )
    if(dtstr_vec=='2')
        V = fread(fid,2*prod(sz),'*single');
        V = reshape(V,[2 sz]);
        V = permute(V,[2 3 1]);
    else
        % read data
        V = fread(fid,prod(sz),'*single');
        V = reshape(V,sz);
    end
  elseif( strcmp(dtstr,'double') > 0 )
    % read data
    V = fread(fid,prod(sz),'*double');
    V = reshape(V,sz);
  end
  
elseif( strcmp(svstr,'VECTORS') > 0 )
  if( strcmp(dtstr,'float') > 0 ) 
    % read data
    V = fread(fid,3*prod(sz),'*single');
    V = reshape(V,[3 sz]);
    V = permute(V,[2 3 4 1]);
  elseif( strcmp(dtstr,'double') > 0 ) 
          V = fread(fid,3*prod(sz),'*double');
    V = reshape(V,[3 sz]);
    V = permute(V,[2 3 4 1]);
  end
end

fclose(fid);