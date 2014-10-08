function [template T subj] = cpd_affine_groupwise(subj, template, T, opt, sigma2)
% T is the initial guess for aligned template to subjects

setNr = length(subj);
D = size(subj{1}, 2);

M = size(template, 1);
N = zeros(length(subj), 1); %Number of surface points per subject
B = cell(1,setNr); %Affine matrix for each member of population
Tr = cell(1,setNr); %Translation vector for each member of population
q = zeros(length(subj), 1); %Objective function (Q) for each member of population 
dQt = zeros(length(subj), 1); %Difference for (Q) that procs template update
dQg = zeros(length(subj), 1); %Difference for (Q) that terminates groupwise registration
sigma2 = zeros(length(subj), 1); %Placeholder for sigma2 for each member of population

for i=1:setNr
    N(i) = size(subj{i}, 1);
end

for i=1:setNr
    sigma2(i)=(M*trace(subj{i}'*subj{i})+N(i)*trace(T{i}'*T{i})-2*sum(subj{i})*sum(T{i})')/(M*N(i)*D);
end

iter=0;

while ( iter < opt.maxGiters && sum(dQg) > opt.QgTol )
    for setInd=1:setNr
        %E-step
        [P, P1, Pt1, Np] = cpd_P(subj{setInd}, T{setInd}, sigma2(setInd), opt.w);
        
        mux = sum(P*subj{setInd},1)/Np;
        muy = sum(P'*T{setInd},1)/Np;
        
        XX = bsxfun(@minus, subj{setInd}, mux);
        YY = bsxfun(@minus, T{setInd}, muy);
        A = XX'*P'*YY;
        
        YPY = (YY'*diag(P1)*YY);
        
        % account for possible singular matrix, checking
        % condition number.  If poorly conditioned,
        % use pseudo-inverse
        if (cond(YPY) < 1e10)
            B{setInd} = A/YPY;
        else
            B{setInd} = A*pinv(YPY);
        end
        Tr{setInd} = mux' - B*(muy');
        
        % output
        T{setInd} = bsxfun(@plus, Y*(B{setInd}'), Tr{setInd}');
        
        qtprev = q;
        
        % estimate error (fast)
        trAB = trace(A*B{setInd}');
        xPx = (Pt1')*sum(XX.*XX,2);
        trBYPYB = trace(B{setInd}*YPY*B{setInd}');
        
        q(setInd) = (xPx-2*trAB+trBYPYB)/(2*sigma2(setInd))+D*Np/2*log(sigma2(setInd));
        dQt(setInd) = abs(q-qtprev);
        
        sigma2(setInd) = (xPx-trAB)/(Np*D);
        
        if (sigma2(setInd) <= 0)
            sigma2(setInd) = opt.QtTol/10; % fix for near zero sigma
        end
        
    end
    
    if ( sum(dQt) < opt.QtTol || iter == (opt.maxGiters-1) )
        
        template = cpd_affine_mean(subj, template, T, A, Tr, sigma2, opt.w);
        
        dQg(setInd) = abs(q-qgprev);
        qgprev = q;
        
        % Update the template
        for setInd=1:setNr
            T{setInd}=template*B{setInd}'+repmat(Tr{setInd}',[M 1]);
        end
    end
   
    iter=iter+1;
    
end

for setInd=1:setNr
    [rot tr] = umeyama(T{setInd}, template);
    subj{setInd} = bsxfun(@plus, rot*subj{setInd}', tr)';  
    T{setInd} = bsxfun(@plus, rot*T{setInd}', tr)'; 
end



