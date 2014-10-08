function [SSM correspondences] = GMM_groupwiseRegistration(subj, opt)
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
for i=1:subjNr
%     tmp = reducepatch(subj{i},opt.pointNr*1.3*2);
%     subjPoint{i} = tmp.vertices;
    subjPoint{i} = refinePoints(subj{i}.vertices,'pointNr', opt.pointNr*1.3);
end

% tmp = reducepatch(subj{floor(subjNr/2)},opt.pointNr*2);
% opt.meanShape = tmp.vertices;
opt.meanShape = refinePoints(subj{i}.vertices,'pointNr', opt.pointNr);

[meanShape correspondences] = groupwiseRegistration(subjPoint, opt);
SSM = generateAtlas(meanShape, correspondences);