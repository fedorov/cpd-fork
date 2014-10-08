function [template T subj] = cpd_fem_groupwise(subj, template, T, opt)

setNr = length(subj);
D = size(subj{1}, 2);
M = size(template, 1);
N = zeros(length(subj), 1); %Number of surface points per subject
qtprev = zeros(length(subj), 1); %Previous objective function (Q) for each member of population 
q = zeros(length(subj), 1); %Objective function (Q) for each member of population 
dQt = zeros(length(subj), 1); %Difference for (Q) that procs template update
dQg = zeros(length(subj), 1); %Difference for (Q) that terminates groupwise registration
sigma2 = zeros(length(subj), 1); %Placeholder for sigma2 for each member of population
%V = cell(1,setNr); %Placeholder for deformation on surface for each member of population
variation = cell(1,setNr);

% FEM stuff
D_material = fem_material_linear(opt.E, opt.nu); 

for i=1:setNr
    N(i) = size(subj{i}, 1);
end

for i=1:setNr
    sigma2(i)=(M*trace(subj{i}'*subj{i})+N(i)*trace(template'*template)-2*sum(subj{i})*sum(template)')/(M*N(i)*D)/10;
end

iter=0;
w=opt.w;
beta = opt.beta;

%magic numbers
mv=0.2; b=2; b_decay=0.7; beta_growth = 10; sigma2_growth = 10;

while ( iter<opt.maxGiters && sum(dQg) > opt.QgTol )

    avg_beta = iter * beta_growth;
    
     parfor setInd=1:setNr
        %E-step
        [P, P1, Pt1, Np] = cpd_P(subj{setInd}, T{setInd}, sigma2(setInd), w);
        
        % M-step
        y_vec = reshape(T{setInd}',[],1);
        [nodes, ~, elems] = tetgen_mex(T{setInd}', F',[],'');
        Phi = [speye(M),zeros(M, size(nodes,2)-M)];
        Phi_tilde = kron(Phi, eye(D));
        fem = fem_model(nodes', elems');
        [K, ~] = getStiffnessMatrix(fem, D_material);
        
        dP1 = spdiags(kron(P1,ones(D,1)),0,D*M,D*M);
        LHS = (Phi_tilde')*dP1*Phi_tilde + beta*sigma2(setInd)*K;
        RHS = -(P*subj{setInd});
        RHS = -Phi_tilde'*(reshape(RHS',[],1)+dP1*(y_vec));
        
        v_vec = LHS\RHS;
        v_vec = Phi_tilde*v_vec;
        V = reshape(v_vec,D,[])';
        
        T{setInd} = template + V;
        variation{setInd} = V{setInd};
        
        sigma2(setInd)=abs((sum(sum(subj{setInd}.^2.*repmat(Pt1,1,D)))+...
            sum(sum(T{setInd}.^2.*repmat(P1,1,D))) ...
            -2*trace(PX'*T{setInd})) /(Np*D));
        
        qtprev(setInd) = q(setInd);
        q(setInd) = cpd_Q(subj{setInd}, T{setInd}, P, P1, Pt1, sigma2(setInd));
        dQt(setInd) = abs(q(setInd)-qtprev(setInd));
        
        if (sigma2 <= 0)
            sigma2 = opt.QtTol/10; % fix for near zero sigma
        end
    end
    
    if(sum(dQt) < ssmOpt.QtTol || iter == (opt.maxGiters-1))
        
        for i=1:setNr
            m_sum = m_sum + T{i};
        end
        template = m_sum / setNr;
        
        for i=1:setNr
            variation{i} = T{i} - template;
        end
        
        r = mean(sigma2)^0.5/2;
        b = b * b_decay;
        disp('Optimize mean shape');
        [template_new T] = cpd_nonrigid_mean(subj, template, variation, sigma2, mv, r, b, avg_beta);
        mv = max(sum((template-template_new).^2, 2).^0.5);
        template = template_new;
        
        sigma2 = sigma2 * sigma2_growth;
    end
    iter=iter+1;
end  