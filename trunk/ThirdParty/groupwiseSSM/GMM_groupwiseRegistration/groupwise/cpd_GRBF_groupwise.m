function [template T subj] = cpd_GRBF_groupwise(subj, template, T, opt, sigma2)

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
for i=1:setNr
    W{i} = zeros(M, D);
end
if(sigma2 == 0)
    for i=1:setNr      
        sigma2(i)=(M*trace(subj{i}'*subj{i})+N(i)*trace(template'*template)-2*sum(subj{i})*sum(template)')/(M*N(i)*D)/10;
    end
end
sigma2 = sigma2 * 8;


iter=0; ntol=opt.tol+10; L=double(ones(1, setNr));
TimesToReset = 5;
G = cpd_G(template,template,opt.beta);

sigma_final = mean(m_d)/10;


check_tol = mean(sigma2)/10;
check_step = 0.5;
mv = 0.2;

b = 2;


while iter<opt.max_it && (ntol > opt.tol) 

    avg_lambda = iter*10; 
    
    if(mean(sigma2) < check_tol && iter > 1)
        check_tol = check_tol * check_step;
        
        
        m_sum = 0;
        for i=1:setNr
            m_sum = m_sum + T{i};
        end
        template = m_sum / setNr;
        for i=1:setNr
            variation{i} = T{i} - template;
        end
        
        
        
        
        
        r = mean(sigma2)^0.5/2;
        b = b * 0.7;
        disp('Optimize mean shape');
        [template_new T] = FindOptimumTempPoint(subj, template, variation, sigma2, mv, r, b, avg_lambda);
        mv = max(sum((template-template_new).^2, 2).^0.5);
        template = template_new;
        
        
        A = template;
        sigma2 = sigma2 * 10;
        G = cpd_G(template,template,opt.beta);
    end
   

    parfor setInd=1:setNr
        L_old(setInd)=L(setInd);

        [P1,Pt1, PX, L(setInd)]=...
            cpd_P(subj{setInd}, T{setInd}, sigma2(setInd) ,0);
        L(setInd)=L(setInd)+opt.lambda/2*trace(W{setInd}'*G*W{setInd});
        ntol=abs((L(setInd)-L_old(setInd))/L(setInd));
        % M-step. Solve linear system for W.

        dP=spdiags(P1,0,M,M); % precompute diag(P)
        W{setInd}=(dP*G+opt.lambda*sigma2(setInd)*eye(M))\(PX-dP*template);
        test = isnan(W{setInd});
        if(test(1,1) == 0)
            T{setInd} = template + G*W{setInd};
            variation{setInd} = G*W{setInd};
        end
        Np=sum(P1);
        
        sigma2(setInd)=abs((sum(sum(subj{setInd}.^2.*repmat(Pt1,1,D)))+...
            sum(sum(T{setInd}.^2.*repmat(P1,1,D))) ...
            -2*trace(PX'*T{setInd})) /(Np*D));    
        if(isnan(sigma2(setInd)))
            sigma2(setInd) = sigma_final;
        end
    end
    disp(['GMM groupwise nonrigid : dL= ' num2str(ntol) ', iter= ' ...
        num2str(iter) ', subject= ' num2str(setNr) ...
        ' sigma2= ' num2str(mean(sigma2))]);


    
    if(opt.viz)
        cpd_plot_iter(template, template);
    end
    iter=iter+1;
    

end
for i=1:setNr
    C{i}=cpd_Pcorrespondence(subj{i},T{i},sigma2(i),0);
end    


