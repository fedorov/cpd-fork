function [T_new T_moved] = cpd_nonrigid_mean(X, T, V, sigma, movingWidth, r, b, lambda)
% I am not even going to try to understand what he did ATM!!!! (Siavash)

[M D] = size(T);
point_nr = 20;

[X_sampled_ind Pr_sampled] = generatePr(X, T, V, sigma, point_nr);
X_all = [];
options = optimset('GradObj','off', 'Display', 'off', 'LargeScale', 'off', 'TolFun', 1.e-3, 'TolX', 1.e-3); %, 'MaxIter', 0);

coef = distDeformedTempTarget(X, T, V, r, b);

m_min = min(T);
m_max = max(T);

Width = max(m_max-m_min);
ctrl_nr = 25;
grindNr_max = 30;
gridNr = min(round(Width/(2*movingWidth)*ctrl_nr^(1/3)), grindNr_max);
ctrl_nr = min(ctrl_nr, gridNr^3);
k_i=2*ones(ctrl_nr,1);
subjNr = length(V);
step = (m_max-m_min)/(gridNr - 1);
if(size(T, 2) == 3)
    [ctrlX ctrlY ctrlZ] = meshgrid((m_min(1):step(1):m_max(1)), ...
        (m_min(2):step(2):m_max(2)), ...
        (m_min(3):step(3):m_max(3)));
    ctrlPt = [ctrlX(:) ctrlY(:) ctrlZ(:)];
    [boundX boundY boundZ] = meshgrid([m_min(1)-1 m_max(1)+1], ...
        [m_min(2)-1 m_max(2)+1], ...
        [m_min(3)-1 m_max(3)+1]);
    bound = [boundX(:) boundY(:) boundZ(:)];
else
    [ctrlX ctrlY] = meshgrid((m_min(1):step(1):m_max(1)), ...
        (m_min(2):step(2):m_max(2)));
    ctrlPt = [ctrlX(:) ctrlY(:)];
    [boundX boundY] = meshgrid([m_min(1)-1 m_max(1)+1], ...
        [m_min(2)-1 m_max(2)+1]);
    bound = [boundX(:) boundY(:)];
end
bound = [];
train_size = 500;
sample_step = max(round(size(T, 1)/train_size), 1);
T_sampled = T(1:sample_step:end,:);
k = 2*ones(floor(size(T_sampled, 1))+size(bound, 1), 1);
for m=1:length(V)
    p_t = [T_sampled; bound];
    v_t = [V{m}(1:sample_step:end,:); zeros(size(bound))];
    [W{m}]=train_rbf(p_t, v_t, p_t, k, 'phs');
    ctrlMv{m} = sim_rbf(p_t, ctrlPt, W{m}, k, 'phs');
end
for i=1:M
    t = T(i,:);
    dist = sum(bsxfun(@minus, ctrlPt, t).^2, 2);
    [tmp ind] = sort(dist);
    SS(i,1) = tmp(ctrl_nr);
    ctrl_sampled = ctrlPt(ind(1:ctrl_nr),:);
    for m=1:subjNr
        V_sampled = ctrlMv{m}(ind(1:ctrl_nr),:);
        W{m}=train_rbf(ctrl_sampled,V_sampled,ctrl_sampled,k_i,'phs');
        X_send{m} = X{m}(X_sampled_ind{m}(i,:), :);
        Pr_send{m} = Pr_sampled{m}(i, :)';
    end
    [T_new(i,:) p] = fminunc(@groupwise_cost, t, options, Pr_send, X_send, ctrl_sampled, W, 2, coef(i,:), lambda);

    
    for m=1:subjNr
        [vt phi]=sim_rbf(ctrl_sampled,T_new(i,:),W{m},k_i,'phs');
        T_moved{m}(i,:) = T_new(i,:)+vt;
    end
end
dist = sum((T_new - T).^2,2);
max((dist ./ SS).^0.5);








function coef = distDeformedTempTarget(X, T, V, r, b)
for i=1:length(V)
    Y = T+V{i};
    tree = kdtree_build( X{i} );
    index = kdtree_nearest_neighbor(tree, Y);
    allD(:,i) = sum((Y-X{i}(index,:)).^2, 2).^0.5;
    kdtree_delete(tree);
%     [tmp, tmp, root] = kdtree(X{i}, []);
%     [index, allD(:,i), root] = kdtreeidx([], Y, root);
%     kdtree([],[],root);
    
end

coef = zeros(size(T, 1), length(V));
for i=1:size(T, 1)
    t = sum(bsxfun(@minus, T, T(i,:)).^2, 2).^0.5;
    d = allD(t < b, :);
    coef(i,:) = exp(-mean(d)/r);
end