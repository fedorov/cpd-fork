function newTemplate = cpd_affine_mean(subj, template, T, A, Tr, sigma2, w)

newTemplate = zeros(size(template));

PX = zeros(size(template));

P1 = cell(1,length(A));

for k=1:length(A)
    movedSubj = A{k}' * bsxfun(@minus, subj{k}, Tr{k}')';
    [Pr, ~, ~, ~] = cpd_P(subj{k}, T{k}, sigma2(k), w);
    P1{k} = sum(Pr, 2);
    PX = PX + Pr*(movedSubj');
end

for j=1:size(template, 1)
    Gamma = zeros(size(A{1}));
    for k=1:length(A)
        Gamma = Gamma + P1{k}(j)*(A{k}'*A{k});
    end
    newTemplate(j,:) = (Gamma)\(PX(j,:)');
end