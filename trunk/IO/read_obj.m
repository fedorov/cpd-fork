function [V, F] = read_obj( filename )
%READ_OBJ Reads a simple OBJ files with v and f only, no relative indexing
%   [V, F] = read_obj(filename)
%   V: Nx3 vertices
%   F: MxK faces
%   filename: input file name

fid = fopen(filename, 'r');

V = [];
F = [];
while (~feof(fid))
   line = nextLine(fid);
   if (line(1) == 'v')
       data = sscanf(line(2:end), '%f', 3);
       V = [V; data'];
   elseif (line(1) == 'f')
       data = sscanf(line(2:end), '%f');
       F = [F; data'];
   end
end

fclose(fid);

    function l = nextLine (fid)
        l = fgetl(fid);
        while (~feof(fid) && (isempty(l) || l(1) == '#'))
            l = fgetl(fid);
        end
    end

end