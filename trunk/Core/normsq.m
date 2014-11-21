function [ out ] = normsquared( in, k )
%NORMSQUARED computes the squared 2-norm along dimension k

    out = sum(in.*in, k);

end

