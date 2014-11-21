function [ fem ] = createBoundingFem( pnts, margins, res )
%CREATEBOUNDINGFEM creates a hex fem beam that bounds the given points
%   fem = createBoundingFem(pnts, margins)
%           creates a beam that bounds the supplied bounds, with extra
%           space around the best-fit axes according to the supplied
%           margins
%   pnts:    Nx3 set of 3D points
%   margins: either 1x1, 3x1, or 6x1 set of margin specifications
%            extra space around best-fit x,y,z axes
%            1x1: same value all around
%            1x3: x, y, z margin
%            1x6: -x,+x,-y,+y,-z,+z margin
%   res:     1x3 number of elements along each dimension
[R, c, w] = tight_box(pnts);

if (nargin < 2 || isempty(margins))
    margins = 0;
end

% correct for margin
mm = zeros(6,1);
if (length(margins) == 1)
    mm = ones(6,1)*margins;
elseif (length(margins) == 3)
    mm = [margins(1), margins(1), margins(2), margins(2), margins(3), margins(3)];
elseif (length(margins) ==6)
    mm = margins;
end

w(1) = w(1)+mm(1)+mm(2);
w(2) = w(2)+mm(3)+mm(4);
w(3) = w(3)+mm(5)+mm(6);

c(1) = c(1) + mm(2)-mm(1);
c(2) = c(2) + mm(4)-mm(3);
c(3) = c(3) + mm(6)-mm(5);

fem = fem_model.createBeam(w, res, c, R);

end

