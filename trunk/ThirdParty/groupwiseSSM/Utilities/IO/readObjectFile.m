function fv = readObjectFile( fileName )
%READOBJECTFILE Summary of this function goes here
%   Detailed explanation goes here
fid = fopen(fileName);
fv.faces = [];
fv.vertices = [];
while ~feof(fid)
    str = fgetl(fid);    
    if(length(str) > 1)
        if(str(1) == 'f')
            fv.faces = [fv.faces; str2num(str(2:end))];
        elseif(str(1) == 'v')
            fv.vertices = [fv.vertices; str2num(str(2:end))];
        end
    end

end
fclose(fid);
end

