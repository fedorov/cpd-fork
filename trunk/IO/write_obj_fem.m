function [] = write_obj_fem( V, E, F, filename, dfmt )
%WRITE_OBJ_FEM Simple routine for writing an OBJ file
%   write_obj(V, E, F, filename, dfmt);
%   V: Nx3  list of vertices
%   E: cell array of elements
%   F: cell array of faces
%   filename: output file name
%   dfmt: format for printing doubles in fprintf (defaults to '% .8f')

if (nargin < 5 || isempty(dfmt))
    dfmt = '% .8f';
end
dfmts = [dfmt ' ']; % with trailing space

fout = fopen(filename,'w');

N = size(V,1);
J = length(E);
K = length(F);


for i=1:N
    fprintf(fout,['v ', dfmts, dfmts, dfmt,'\n'], V(i,1), V(i,2), V(i,3));
end

for i=1:J
    fprintf(fout,'e ');
    elem = E{i};
    n = size(elem,2);
    if (n==8)
        elem = elem([4:-1:1,8:-1:5]);
    end
    for j=1:(n-1)
        fprintf(fout,'%i ', elem(j));
    end
    fprintf(fout,'%i\n', elem(end));
end

for i=1:K
    fprintf(fout,'f ');
    f = F{i};
    n = size(f,2);
    for j=1:(n-1)
        fprintf(fout,'%i ', f(j));
    end
    fprintf(fout,'%i\n', f(end));
end

fclose(fout);