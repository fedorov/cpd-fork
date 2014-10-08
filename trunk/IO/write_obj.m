function [] = write_obj( V, F, filename, dfmt )
%WRITE_OBJ Simple routine for writing an OBJ file
%   write_obj(V, F, filename, dfmt);
%   V: Nx3  list of vertices
%   F: MxK  triangular faces
%   filename: output file name
%   dfmt: format for printing doubles in fprintf (defaults to '% .8f')

if (nargin < 4 || isempty(dfmt))
    dfmt = '% .8f';
end
dfmts = [dfmt ' ']; % with trailing space

fout = fopen(filename,'w');

N = size(V,1);
M = size(F, 1);
K = size(F, 2);

for i=1:N
    fprintf(fout,['v ', dfmts, dfmts, dfmt,'\n'], V(i,1), V(i,2), V(i,3));
end

for i=1:M
    fprintf(fout,'f ');
    for j=1:(K-1)
        fprintf(fout,'%i ', F(i,j));
    end
    fprintf(fout,'%i\n', F(i,end));
end

fclose(fout);