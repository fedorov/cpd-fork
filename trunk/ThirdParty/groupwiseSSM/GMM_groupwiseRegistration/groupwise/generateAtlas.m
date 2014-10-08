function [atlas deformed] = generateAtlas(template, deformed)
setNr = length(deformed);
[M D] = size(template);
obj_vari = zeros(setNr, D*M);
for i=1:setNr
    for setInd=1:setNr
        [rot tr] = umeyama(double(deformed{setInd}), double(template));
        deformed{setInd} = bsxfun(@plus, rot*double(deformed{setInd}'), tr)';  
    end    
    tmp = double(deformed{i}') - double(template');
    obj_vari(i, :) = tmp(:)';
end
[coef, score, latent] = princomp(obj_vari, 'econ');
atlas.mod_nr = length(latent);
meanVar = mean(obj_vari);
atlas.mean = template+reshape(meanVar', D, length(meanVar)/D)';
atlas.mean = bsxfun(@minus, atlas.mean, mean(atlas.mean));
atlas.latent = latent;
atlas.mods = coef;
atlas.score = score;