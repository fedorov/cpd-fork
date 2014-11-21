function [meanShape var] = cpd_groupwise(subj, opt)
%CPD_GROUPWISE performs the groupwise registration and mean shape
%calculation
%   Input:
%   subj    1xN or Nx1 cell vector of population
%   opt     Options for SSM construction
%
%   Output:
%   meanShape   Mx3 array of mean shape surface points
%   var         1xN cell of population variations w.r.t the mean shape


setNr = length(subj);
meanShape = opt.meanShape;
meanShape = bsxfun(@minus, meanShape, mean(meanShape));
meanShape = bsxfun(@rdivide, meanShape, std(meanShape));
subj_std = cell(1,setNr);
var = cell(1,setNr);

% Rigid registration of population to mean (initial alignment)
for i=1:setNr
    subj{i} = bsxfun(@minus, subj{i}, mean(subj{i}));
    subj_std{i} = std(subj{i});
    subj{i} = bsxfun(@rdivide, subj{i}, subj_std{i});
    [ TY, ~, ~, ~, ~, ~ ] = cpd_rigid( meanShape, subj{i}, opt.w, opt.tol, opt.maxiters, [], [], [], []);
    subj{i} = TY;
    i
end

for i=1:setNr
    var{i} = meanShape;
end
disp('Template initialized');

% Affine registration of the mean to the population
if(opt.affine == 1)
    [~, var, subj] = cpd_affine_groupwise(subj, meanShape, var, opt, 0);
    S = 0;
    for i=1:setNr
        S = S + var{i};
    end
    meanShape = S./setNr;
    if(opt.save == 1)
        save('affine.mat', 'meanShape', 'var', 'subj');
    end
end

if(opt.save == 1)
    save('meanShapeAndVarEMICP.mat', 'meanShape', 'var');
end

% Non-rigid registration of the mean to the population
if(opt.nonrigid == 1)
    [~, var, ~] = cpd_fem_groupwise(subj, meanShape, var, opt);
    S = 0;
    for i=1:setNr
        S = S + var{i};
    end
    meanShape = S./setNr;
end

for i=1:setNr
    var{i} = bsxfun(@times, var{i}, subj_std{i});
end


if(opt.save == 1)
    save('meanShapeAndVarGRBF.mat', 'meanShape', 'var');
end


end

