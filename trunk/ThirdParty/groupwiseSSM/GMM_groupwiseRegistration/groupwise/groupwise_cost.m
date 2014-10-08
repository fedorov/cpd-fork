function c = groupwise_cost(t, pr, X, T, W, s, coef, lambda)
[M D] = size(T);
k_i=s*ones(M,1);
c = 0;
% lambda = 50;%20;
vt_mean = zeros(size(t));
for m=1:length(pr)
    if(s == 3)
        [vt phi]=sim_rbf(T,t,W{m},k_i,'gaussian');
    else
        [vt phi]=sim_rbf(T,t,W{m},k_i,'phs');
    end
    vt_mean = vt_mean + vt;
    new_t = t+vt;
    dist = sum(bsxfun(@minus, X{m}, new_t).^2,2).^0.5;
    c = c + sum(pr{m}.*dist)*coef(m);
end
vt_mean = vt_mean / length(pr);
c = c + sum(vt_mean.^2)*lambda;
