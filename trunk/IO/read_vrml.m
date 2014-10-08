function [vertex, face, vert_norm, face_norm] = read_vrml(filename)

% read_wrl - load a mesh from a VRML file
%
%   [vertex, face] = read_wrl(filename);
%
%   Copyright (c) 2004 Gabriel Peyr

fp = fopen(filename,'r');
if fp == -1
    fclose all;
    error(['Cannot open file ' filename]);
end

tempstr = ' ';

key = 'point [';
vertex = [];

while ( tempstr ~= -1)
    tempstr = fgets(fp);   % -1 if eof
    if( ~isempty(findstr(tempstr,key)) )
        [vertex,nc] = fscanf(fp,'%f%c %f%c %f%c', Inf);
        
        vertex = reshape(vertex, [6,length(vertex)/6]);
        vertex([2 4 6],:) = [];
        
        if 0 %% old code
            nc = 3;
            while nc>0
                tempstr = fgets(fp);   % -1 if eof
                [cvals,nc] = sscanf(tempstr,'%f %f %f,');
                if nc>0
                    if mod(nc,3)~=0
                        error('Not correct WRL format');
                    end
                    vertex = [vertex, reshape(cvals, 3, nc/3)];
                end
            end
        end
        
        break;
    end
end

key = 'vector [';
vert_norm = [];

while (( tempstr ~= -1) & ( nargout > 2 ))
    tempstr = fgets(fp);   % -1 if eof
    if( ~isempty(findstr(tempstr,key)) )
        [vert_norm,nc] = fscanf(fp,'%f%c %f%c %f%c', Inf);
        
        vert_norm = reshape(vert_norm, [6,length(vert_norm)/6]);
        vert_norm([2 4 6],:) = [];
        
        if 0 %% old code
            nc = 3;
            while nc>0
                tempstr = fgets(fp);   % -1 if eof
                [cvals,nc] = sscanf(tempstr,'%f %f %f,');
                if nc>0
                    if mod(nc,3)~=0
                        error('Not correct WRL format');
                    end
                    vert_norm = [vert_norm, reshape(cvals, 3, nc/3)];
                end
            end
        end
        
        break;
    end
end

key = 'coordIndex [';
face = [];

while (( tempstr ~= -1) & ( nargout > 1))
    tempstr = fgets(fp);   % -1 if eof
    if( ~isempty(findstr(tempstr,key)) )
        
        % should work with , and space separators
        [face,nc] = fscanf(fp,'%d%c %d%c %d%c -1%c', Inf);
        face = reshape(face, [7,length(face)/7])+1;
        face([2 4 6 7],:) = [];
        
        if 0 %% old code
            nc = 3;
            while nc>0
                tempstr = fgets(fp);   % -1 if eof
                [cvals,nc] = sscanf(tempstr,'%d %d %d -1,');
                if nc==0 || mod(nc,3)~=0
                    [cvals,nc] = sscanf(tempstr,'%d, %d, %d, -1,');
                end
                if nc>0
                    face = [face, reshape(cvals, 3, nc/3)+1];
                end
            end
        end
        
        fclose(fp);
        tr = TriRep(face', vertex');
        face_norm = faceNormals(tr)';
        
        return;
    end
end
