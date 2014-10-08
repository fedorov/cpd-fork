function [V, E, F] = read_obj_fem( filename )
%READ_OBJ Reads a simple OBJ files with {v, e, f} only, no relative indexing
%   [V, E, F] = read_obj_fem(filename)
%   V: Nx3, vertices
%   E: elements as cell array
%   F: faces as cell array
%   filename: input file name

fid = fopen(filename, 'r');

V = [];
E = {};
F = {};
while (~feof(fid))
   line = nextLine(fid);
   if (line(1) == 'v')
       data = sscanf(line(2:end), '%f', 3);
       V = [V; data'];
   elseif (line(1) == 'e')
       data = sscanf(line(2:end), '%f');
       % reverse hex due to ArtiSynth inconsistency
       if (numel(data)==8)
           data = data([4:-1:1,8:-1:5]);
       end
       E = [E; {data'}];
   elseif (line(1) == 'f')
       data = sscanf(line(2:end), '%f');
       F = [F; {data'}];
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