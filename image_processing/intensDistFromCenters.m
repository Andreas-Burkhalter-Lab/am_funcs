%%% for each pix in projim, get distance to nearest center in
%%% patchcentersim
% results = intensDistFromCenters(patchcentersim,projim,zoom,scope)

function intensdist = intensDistFromCenters(patchcentersim,projim,roi,zoom,scope)

patchcentersim = loadbw(patchcentersim);
roi = loadbw(roi);

if ischar(projim)
    projim = imread(projim);
end
projim = projim(:,:,1);

if ~exist('scope','var')
    scope = [];
end

if any(size(projim) ~= size(patchcentersim))
    error('image dimension mismatch')
end
if any(size(roi) ~= size(patchcentersim))
    error('image dimension mismatch')
end

patchcentersim_roi = patchcentersim;
patchcentersim_roi(~roi) = false;
[patchpixy patchpixx] = find(patchcentersim); % all pixels in patches

% get centers of mass
[bounds, labelmat, npatches] = bwboundaries(patchcentersim_roi);
patchcentersYX = NaN(npatches,2);
for i = 1:npatches
    [ycenter xcenter] = find(labelmat==i);
    patchcentersYX(i,:) = [mean(ycenter) mean(xcenter)];
end

[y x] = find(roi);
npix = length(x);
intensdist = table(x,y,NaN(npix,1),NaN(npix,1),NaN(npix,1),'VariableNames',{'x','y','closest_neighbor','distcenter_pix','projintens'});
intensdist.projintens = projim(sub2ind(size(projim),y,x));
[intensdist.closest_center intensdist.distcenter_pix] = knnsearch(patchcentersYX,[y x],'K',1); % dist from nearest center
[~, intensdist.distpatch_pix] = knnsearch([patchpixy patchpixx],[y x],'K',1);
intensdist.distcenter_um = intensdist.distcenter_pix * umPerPix(zoom,scope); % dist from nearest patch pixel
intensdist.distpatch_um = intensdist.distpatch_pix * umPerPix(zoom,scope);


