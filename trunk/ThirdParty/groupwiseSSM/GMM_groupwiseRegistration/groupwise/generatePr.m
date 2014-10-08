function [X_sampled_ind Pr_sampled] = generatePr(X, T, V, sigma, point_nr)
% X is a set
if(~exist('tmp','dir'))
    mkdir('tmp');
end
for i=1:length(X)
    m_Pr = cpd_Pr(X{i}, T+V{i}, sigma,0,ones(size(X{i},1),1));
    [A I] = sort(m_Pr, 2, 'descend');
    m_mean = mean(A);
    for j=1:size(A, 2)
        l(j) = sum(m_mean(1:j));
    end
    point_nr = size(A, 2);
    for j=1:size(A, 2)
        if(sum(m_mean(1:j))/sum(m_mean) > 0.9)
            point_nr = j;
            break;
        end
    end    
    for j=1:size(A, 2)
        if(sum(m_mean(1:j))/sum(m_mean) > 0.999)
            point_nr = j;
            break;
        end
    end        
    clear m_pr;
    Pr_sampled{i} = A(:,1:point_nr);
    X_sampled_ind{i} = I(:, 1:point_nr);
end
