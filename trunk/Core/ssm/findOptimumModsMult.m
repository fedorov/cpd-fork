function b = findOptimumModsMult(A, S, Pr, mods, m_latent, eigenCoefWeight)
%   abbrevations:
%   R:      Number of mods
%   A:      Atlas points: MxD
%   S:      Subject points: NxD
%   Pr:     Probability MxN
%   meanVar:Mean variation of atlas MxD
%   mods:   mods of variation DMxR
[M D] = size(A);
[N D] = size(S);

R = size(mods, 2);
Gam = zeros(1, R);
for i=1:M
    f = bsxfun(@minus, S, A(i,:));
    res = f*mods((i-1)*D+1:i*D,:);
    %Gam = Gam + sum(bsxfun(@times, res, Pr(i, :)'));
    Gam = Gam + Pr(i, :)*res;
end

Ups = zeros(R, R);
for i=1:M
    res = mods((i-1)*D+1:i*D,:)' * mods((i-1)*D+1:i*D,:);
    Ups = Ups + sum(Pr(i,:))*res;
end
b = (Gam/(Ups+eigenCoefWeight*diag(1./(m_latent.^0.5))))';


clear Gam Ups
