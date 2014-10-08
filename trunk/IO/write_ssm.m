function [] = write_ssm( X, P, S, F, filename, dfmt )
%WRITE_SSM Simple routine for writing an SSM to a file
%   write_ssm(X, P, S, F, filename);
%   X: NxD mean shape
%   P: D*NxK modes of variation
%   S: Kx1  eigenvalues
%   F: Mx3  triangular faces
%   filename: output file name
%   dfmt: format for printing doubles in fprintf (defaults to '% .8f')

% output file:
% # Points, Dimension, Modes, Faces
% # eigenvalues
% # Mean, modes
% # Faces, indexed starting at 0

if (nargin < 6 || isempty(dfmt))
    dfmt = '% .8f';
end
dfmts = [dfmt ' ']; % with trailing space

fout = fopen(filename,'w');

N = size(X,1);
D = size(X,2);
K = size(P, 2);
M = size(F, 1);

fprintf(fout, '# SSM File Format 0.1\n');

fprintf(fout, '# Points, Dimension, Modes, Faces\n');
fprintf(fout, '%d %d %d %d\n\n', N, D, K, M);

fprintf(fout, '# Eigenvalues\n');

for i=1:K
    fprintf(fout, dfmts, S(i));
end
fprintf(fout,'\n\n');

fprintf(fout, '#Mean, modes\n');
for i=1:(D*N)
    fprintf(fout, dfmts, X(i));
    for j=1:K
        fprintf(fout, dfmts, P(i,j));
    end
    fprintf(fout, '\n');
end
fprintf(fout, '\n\n');

fprintf(fout,'# faces\n');
for i=1:M
    fprintf(fout, '%d %d %d\n', F(i,1)-1, F(i, 2)-1, F(i,3)-1);
end
fprintf(fout, '\n\n');

fclose(fout);

end

