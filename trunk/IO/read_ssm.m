function [ X, P, S, F ] = read_ssm( filename )
%READ_SSM Reads an SSM from a simple file
%   [X, P, S, F] = read_ssm(filename)
%   X: NxD mean shape
%   P: D*NxK modes of variation
%   S: Kx1  eigenvalues
%   F: Mx3  triangular faces
%   filename: intput file name

% input file:
% # Points, Dimension, Modes, Faces
% # eigenvalues
% # Mean, modes
% # Faces, indexed starting at 0

fid = fopen(filename, 'r');

% # Points, Dimension, Modes, Faces
line = nextLine(fid);
data = sscanf(line,'%d', 4);
N = data(1);
D = data(2);
K = data(3);
M = data(4);

% Eigenvalues
line = nextLine(fid);
data = sscanf(line, '%f', K);
S = data(:);

% Mean, modes
X = zeros(N, D);
P = zeros(D*N, K);
for i=1:(D*N)
    line = nextLine(fid);
    data = sscanf(line, '%f', K+1);
    X(i) = data(1);
    P(i,:) = data(2:end);
end

% Faces
F = zeros(M, 3);
for i=1:M
    line = nextLine(fid);
    data = sscanf(line, '%d', 3);
    F(i,:) = data(1:3)+1;
end

fclose(fid);

    function l = nextLine (fid)
        l = fgetl(fid);
        while (~feof(fid) && (isempty(l) || l(1) == '#'))
            l = fgetl(fid);
        end
    end

end

