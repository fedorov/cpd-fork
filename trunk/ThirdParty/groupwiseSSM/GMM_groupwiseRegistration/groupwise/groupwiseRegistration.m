function [meanShape var] = ...    groupwiseRegistration(subj, opt)setNr = length(subj);meanShape = [];% for i=1:setNr%     subj{i} = refinePoints(X{i}, 'pointNr', N);%     meanShape = [meanShape; subj{i}];   % end% meanShape = refinePoints(meanShape, 'pointNr', M);meanShape = opt.meanShape;meanShape = bsxfun(@minus, meanShape, mean(meanShape));meanShape = bsxfun(@rdivide, meanShape, std(meanShape));rigidOpt.method = 'rigid';rigidOpt.scale = 0;rigidOpt.normalize = 0;rigidOpt.tol = 1.e-5;rigidOpt.sigma2 = 100;for i=1:setNr    close all;    subj{i} = bsxfun(@minus, subj{i}, mean(subj{i}));    subj_std{i} = std(subj{i});    subj{i} = bsxfun(@rdivide, subj{i}, subj_std{i});    trans{i} = cpd_register(meanShape, subj{i}, rigidOpt);    subj{i} = trans{i}.Y;end% [IDX meanShape] = kmeans(meanShape, M);% qq = floor(rand*length(X))+1;% meanShape = refinePoints(X{qq}, 'pointNr', M);for i=1:setNr    var{i} = meanShape;enddisp('Template initialized');if(opt.affine == 1)    [meanShape var subj] = cpd_affine_groupwise(subj, meanShape, var, opt, 0);    S = 0;    for i=1:setNr        S = S + var{i};    end    meanShape = S./setNr;    if(opt.save == 1)        save('affine.mat', 'meanShape', 'var', 'subj');    endendif(opt.save == 1)    save('meanShapeAndVarEMICP.mat', 'meanShape', 'var');endif(opt.nonrigid == 1)    [meanShape var subj] = cpd_GRBF_groupwise(subj, meanShape, var, opt, 0);    S = 0;    for i=1:setNr        S = S + var{i};    end    meanShape = S./setNr;endfor i=1:setNr    var{i} = bsxfun(@times, var{i}, subj_std{i});endif(opt.save == 1)    save('meanShapeAndVarGRBF.mat', 'meanShape', 'var');end