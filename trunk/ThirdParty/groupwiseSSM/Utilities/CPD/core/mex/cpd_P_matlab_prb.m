function [P1, Pt1, PX] = cpd_P_matlab_PR( ...
    x, xPr, y, sigma2, outlier)



[M, D] = size(y);
N = size(x, 1);

% PX = zeros(M, 3);
% P1 = zeros(M, 1);
% Pt1 = zeros(N, 1);

% P = zeros(M, 1);


ksig = -2.0*sigma2;
outlier_tmp=(outlier*M*((-ksig*pi)^(0.5*D)))/((1-outlier)*N); 

% 
% for n=1:N
%     P = exp(sum(bsxfun(@minus, y, x(n,:)).^2, 2)/ksig);
%     sp = sum(P)+outlier_tmp;
%     Pt1(n)=1-outlier_tmp/sp*xPr(n);
%     P = P ./ sp;
%     P1 = P1 + P*xPr(n);
%     PX = PX + P*x(n,:)*xPr(n);    
% end
% toc
% 
% 
% 
% 
% tic;
% m_size=norm(max(y)-min(y));
% portion=m_size.^2/size(y,1);
% k_numb = ceil(3*portion*sqrt(sigma2));
% tree = kdtree_build( y );
PX = zeros(M, 3);
P1 = zeros(M, 1);
Pt1 = zeros(N, 1);



k_numb = floor(ones(3,1)*M/10);
for l=1:3
     P = sort(exp(sum(bsxfun(@minus, y, x(floor(rand*N)+1,:)).^2, 2)/ksig), 1, 'descend');
     P_sum = sum(P);
     for k=1:M
         if(sum(P(1:k))/P_sum > 0.99)
             k_numb(l) = k;
             break;
         end
     end

end
k_numb = floor(mean(k_numb));

ns = createns(y,'nsmethod','kdtree');

idx = knnsearch(ns,x,'k',k_numb);

for n=1:N
%     idx = knnsearch(ns,x(n,:),'k',k_numb);
%     idx = kdtree_k_nearest_neighbors(tree,x(n,:),k_numb);
    P = exp(sum(bsxfun(@minus, y(idx(n,:),:), x(n,:)).^2, 2)/ksig);
    sp = sum(P)+outlier_tmp;
    Pt1(n)=(1-outlier_tmp/sp)*xPr(n);
    P = P ./ sp;
    P1(idx(n,:)) = P1(idx(n,:)) + P*xPr(n);
    PX(idx(n,:),:) = PX(idx(n,:),:) + P*x(n,:)*xPr(n);       
end

% kdtree_delete(tree);



% 
% fac_size = 1000;
% 
% for n=1:floor(N/fac_size)
%     xTmp{n} = x((n-1)*fac_size+1:n*fac_size,:);
% end
% 
% parfor n=1:floor(N/fac_size)
% 
%     [P1tmp{n} Pt1tmp{n} PXtmp{n}] = computeP(y, xTmp{n},sigma2, outlier);
% 
% end
% 
% Pt1 = [];
% for n=1:floor(N/fac_size)
%     P1 = P1 + P1tmp{n};
%     PX = PX + PXtmp{n};
%     Pt1 = [Pt1; Pt1tmp{n}];
% end
% % P1 = sum(P1, 2);
% % PX = sum(PX, 3);
% 
% 
% function [P1 Pt1 PX] = computeP(y, x, sigma2, outlier)
%     [M D] = size(y);
% N = size(x, 1);
% ksig = -2.0*sigma2;
% outlier_tmp=(outlier*M*((-ksig*pi)^(0.5*D)))/((1-outlier)*N); 
% PX = zeros(M, 3);
% P1 = zeros(M, 1);
% Pt1 = zeros(N, 1);
% 
% P = zeros(M, 1);
% for n=1:N
%     P = exp(sum(bsxfun(@minus, y, x(n,:)).^2, 2)/ksig);
%     sp = sum(P)+outlier_tmp;
%     Pt1(n)=1-outlier_tmp/sp;
%     P = P ./ sp;
%     P1 = P1 + P;
%     PX = PX + P*x(n,:);    
% end


