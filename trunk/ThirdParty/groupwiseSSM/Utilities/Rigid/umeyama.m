function [R T Z] = umeyama(X, Y)
% UMEYAMA finds the transformation from X to Y
% usage:
%   [R T] = umeyama(X, Y)
Xmean = mean(X);
Ymean = mean(Y);
Xnorm = bsxfun(@minus, X, Xmean);
Ynorm = bsxfun(@minus, Y, Ymean);
[U D V] = svd(Xnorm'*Ynorm);
R = V*U';
T = Ymean' - R*Xmean';

Z = bsxfun(@plus, X*R', T');