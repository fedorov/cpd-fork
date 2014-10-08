function [P1 Pt1 PX E] = cpd_P_matlab( ...
    x, y, sigma2, outlier)

[M D] = size(y);
N = size(x, 1);
ksig = -2.0*sigma2;
outlier_tmp=(outlier*M*((-ksig*pi)^(0.5*D)))/((1-outlier)*N); 



if(~exist('createns', 'file'))
    numb = ones(3,1)*M/10;
    for l=1:3
        P = sort(exp(sum(bsxfun(@minus, y, x(floor(rand*N)+1,:)).^2, 2)/ksig), 1, 'descend');
        P_sum = sum(P);
        for k=1:M
            if(sum(P(1:k))/P_sum > 0.99)
                numb(l) = k;
                break;
            end
        end

    end
    k_numb = ceil(mean(numb));
    P = zeros(M,N);
    ns = createns(y);
    [idx, dist] = knnsearch(ns,x,'k',k_numb);
    idx = bsxfun(@plus, idx',size(P,1)* (0:size(P,2)'-1));
    P(idx) = exp(dist'.^2/ksig);
    E = -sum(log(sum(P)+outlier_tmp));  
    P = bsxfun(@rdivide, P, sum(P)+outlier_tmp);
    Pt1 = sum(P, 1)';
    P1 = sum(P, 2);
    PX = P * x;
    E = E + D*N*log(sigma2)/2;   
else
    PX = zeros(M, 3);
    P1 = zeros(M, 1);
    Pt1 = zeros(N, 1); 
    E = 0;
    for n=1:N

        P = exp(sum(bsxfun(@minus, y, x(n,:)).^2, 2)/ksig);
        sp = sum(P)+outlier_tmp;
        P = P ./ sp;
        Pt1(n)=1-outlier_tmp/sp;
        P1 = P1 + P;
        PX = PX + P*x(n,:);    
        E = E - log(sp);
    end
    E = E + D*N*log(sigma2)/2;  
end