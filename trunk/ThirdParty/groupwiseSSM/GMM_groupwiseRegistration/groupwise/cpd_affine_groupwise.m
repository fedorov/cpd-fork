function [template T subj] = cpd_affine_groupwise(subj, template, T, opt, sigma2)
% T is the initial guess for aligned template to subjects


setNr = length(subj);
D = size(subj{1}, 2);

M = size(template, 1);
N = zeros(length(subj), 1);

m_d = [];
for i=1:setNr
    N(i) = size(subj{i}, 1);
    for l=1:5
        index = floor(rand * size(subj{i}, 1))+1;
        dist = sum(bsxfun(@minus, subj{i}, subj{i}(index, :)).^2, 2);
        [a e] = sort(dist);
        m_d = [m_d; a(2)];
    end
end


if(sigma2 == 0)
    for i=1:setNr
%         sigma2(i)=(M*trace(subj{i}'*subj{i})+N(i)*trace(template'*template)-2*sum(subj{i})*sum(template)')/(M*N(i)*D);
        sigma2(i)=(M*trace(subj{i}'*subj{i})+N(i)*trace(T{i}'*T{i})-2*sum(subj{i})*sum(T{i})')/(M*N(i)*D);
    end    
end
iter=0; ntol=opt.tol+10; L=double(ones(1, setNr));
step = 30;
step_size = 4;

sigma_initial = sigma2;
sigma_factor = 0.95;
sigma_final = 0.005;

while (iter<opt.max_it) && (ntol > opt.tol) && sum(sigma2(1) > 1e-8) > 0
    sigma2 = max(sigma_final, sigma2);
    for setInd=1:setNr
        L_old(setInd)=L(setInd);

        [P1,Pt1, PX, L(setInd)]=cpd_P(subj{setInd}, T{setInd}, sigma2(setInd) ,0); st='';
        ntol=abs((L(setInd)-L_old(setInd))/L(setInd));
        
        Np=sum(P1);
        mu_subj=subj{setInd}'*Pt1/Np;
        mu_temp=template'*P1/Np;

        % Solve for parameters
        B1=PX'*template-Np*(mu_subj*mu_temp');
        B2=(template.*repmat(P1,1,D))'*template-Np*(mu_temp*mu_temp');
        B=B1/B2; % B= B1 * inv(B2);
        t=mu_subj-B*mu_temp;     
        
        A{setInd} = B;
        Tr{setInd} = t;


        sigma2(setInd)=abs(sum(sum(subj{setInd}.^2.*repmat(Pt1,1,D)))- Np*(mu_subj'*mu_subj) -trace(B1*B'))/(Np*D); 
        % Update centroids positioins
        T{setInd}=template*B'+repmat(t',[M 1]);    
        
        variation{setInd} = template - T{setInd};
    end
    step = step - 1;

    if(opt.viz && mod(iter, 10) == 0)
        cpd_plot_iter(template, template);
    end
    
    if(step == 0) 
        step = step_size;
        step_size = ceil(step_size * 1.1);
        
        
        template = FindOptimumTempPointInAffine(subj, template, T, A, Tr, sigma2);
        for setInd=1:setNr
            T{setInd}=template*A{setInd}'+repmat(Tr{setInd}',[M 1]);
        end

    end
   
    disp(['GMM Groupwise affine ' st ' : dL= ' num2str(ntol) ', iter= ' num2str(iter) ', subject= ' num2str(setInd) ' sigma2= ' num2str(mean(sigma2))]);

    if(opt.viz && mod(iter, 1) == 0)
        cpd_plot_iter(template, template);
    end
    iter=iter+1;

    
end
for setInd=1:setNr
%             T{setInd}=template*A{setInd}'+repmat(Tr{setInd}',[M 1]); 
    [rot tr] = umeyama(T{setInd}, template);
    subj{setInd} = bsxfun(@plus, rot*subj{setInd}', tr)';  
%             T{setInd} = template;
    T{setInd} = bsxfun(@plus, rot*T{setInd}', tr)'; 
end
disp('CPD registration succesfully completed.');


