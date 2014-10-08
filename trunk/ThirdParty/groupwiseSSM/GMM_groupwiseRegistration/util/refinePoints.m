function [newP ind Y] = refinePoints(p, varargin)
% parsing the input arguments
parser = inputParser;    
parser.addParamValue('minDist', 0, @(x)isnumeric(x));
parser.addParamValue('pointNr', size(p, 1), @(x)isnumeric(x));
parser.addParamValue('faces', 0, @(x)isnumeric(x));
parser.addParamValue('includeInd', 0, @(x)isnumeric(x));
parser.parse(varargin{:});
minDist = parser.Results.minDist;
pointNr = parser.Results.pointNr;
fv.faces = parser.Results.faces;
fv.vertices = p;
ind = parser.Results.includeInd;

if(fv.faces ~= 0)
    b = reducePatch(fv, pointNr / 2);

    neighbours = zeros(size(p, 1), 10);
    index = ones(size(p, 1), 1);
    for i=1:size(fv.faces, 1)
        neighbours(fv.faces(i, 1), index(fv.faces(i, 1))) = fv.faces(i, 2);
        neighbours(fv.faces(i, 1), index(fv.faces(i, 1))+1) = fv.faces(i, 3);
        index(fv.faces(i, 1)) = index(fv.faces(i, 1)) + 2;
        neighbours(fv.faces(i, 2), index(fv.faces(i, 2))) = fv.faces(i, 1);
        neighbours(fv.faces(i, 2), index(fv.faces(i, 2))+1) = fv.faces(i, 3);
        index(fv.faces(i, 2)) = index(fv.faces(i, 2)) + 2;
        neighbours(fv.faces(i, 3), index(fv.faces(i, 3))) = fv.faces(i, 1);
        neighbours(fv.faces(i, 3), index(fv.faces(i, 3))+1) = fv.faces(i, 2);
        index(fv.faces(i, 3)) = index(fv.faces(i, 3)) + 2;
    end
    for i=1:size(neighbours, 1)
        q = unique(neighbours(i, 1:index(i)-1));
        neighbours(i, :) = 0 * neighbours(i, :);
        index(i) = length(q);
        neighbours(i, 1:index(i)) = q;
    end
    distance = ones(size(p, 1), 1)*1.e10;
    ind = find(sum(bsxfun(@minus, fv.vertices,b.vertices(1,:)), 2) == 0);
    newP = b.vertices(1,:);
    while(size(newP, 1) < pointNr)
        queue = ind(end);
        distance(ind(end)) = 0;
        while(~isempty(queue))
            check = queue(1);
            for m=1:index(check)
                newDist = distance(check)+...
                    norm(p(check,:)-p(neighbours(check, m),:));
                if(newDist < distance(neighbours(check, m)))
                    distance(neighbours(check, m)) = newDist;
                    queue = [queue neighbours(check, m)];
                end
            end
            queue(1) = [];
        end
        if(size(newP, 1) < size(b.vertices,1))
            ind = [ind; find(sum(bsxfun(@minus, fv.vertices,b.vertices(size(newP, 1)+1,:)), 2) == 0)];
            newP = [newP; b.vertices(size(newP, 1)+1,:)];
        else
            tmp_distance = distance .* (distance < 1.e10);
            [Y,I] = max(tmp_distance);
            ind = [ind; I];
            newP = [newP; p(I, :)];        
        end
    end
else
    if(ind == 0)
        ind = floor(rand*size(p,1))+1;
    end
    newP = p(ind,:);
    dist = 1000*ones(size(p,1),1);
    for i=1:length(ind)
        dist = min([dist sum(bsxfun(@minus, p, newP(i,:)).^2, 2).^0.5]')';
    end
    Y = Inf;
    while(Y > minDist && size(newP, 1) < pointNr)
        if(mod(size(newP, 1), 2000) == 0)
            fprintf('\nnew size=%d, min distance=%f', size(newP, 1), Y);
        end
        [Y,I] = max(dist);
        ind = [ind; I];
        newP = [newP; p(I, :)];
        newDist = sum(bsxfun(@minus, p, p(I,:)).^2, 2).^0.5;
        dist = min([dist newDist]')';
        if(size(newP, 1) == size(p, 1))
            break;
        end
    end   
end
