function [SSM correspondences] = cpd_ssmCon(subj, opt)
% subj is a cell which each element is a struct. The struct has two
% elements, vertices and faces. 
% opt is the options for the registration
%   affine      to perform groupwise affine registration
%   nonrigid    to perform groupwise non-rigid registration
%   lambda      non-rigid registration parameter. greater lambda is the
%               more constrained deformation
%   beta        non-rigid registration parameter. greater beta is the
%               more constrained deformation
%   tol         tolerance to stop the registration
%   max_it      maximum iteration of the registration
%   viz         vizualization of the mean shape
%   pointNr     number of the point in the mean shape (max=1500)

subjNr = length(subj);
subjPoint = cell(1,subjNr);

for i=1:subjNr
    subjPoint{i} = refinePoints(subj{i}.vertices,'pointNr', opt.pointNr*1.3);
end

opt.meanShape = refinePoints(subj{i}.vertices,'pointNr', opt.pointNr);

[meanShape correspondences] = cpd_groupwise(subjPoint, opt);
SSM = generateAtlas(meanShape, correspondences);

end

