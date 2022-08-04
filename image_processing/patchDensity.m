 %%% get center-to-nearest-center spacing of w points in a bw image
 % last updated 19/12/15 on thermaltake
 
 function [analysis] = patchDensity(patchesFile,roifile,zoom,scope)
 
 plotting = 0;
 
 densrecov_intervals_um = 10:10:350; % bins for distances for density recovery profile
 
[patchimage analysis.patchesFile] = loadbw(patchesFile);
[roi analysis.roifile] = loadbw(roifile);
if ~exist('scope','var')
    scope = [];
end

if any(size(patchimage) - size(roi) ~=0)
    error('image dimension mismatch')
end
patchesimage_roi = patchimage;
patchesimage_roi(~roi) = false;

% get centers of mass
[bounds, labelmat, npatches] = bwboundaries(patchesimage_roi);
patchcentersYX = NaN(npatches,2);
for i = 1:npatches
    [y x] = find(labelmat==i);
    patchcentersYX(i,:) = [mean(y) mean(x)];
end
    
% nearest neighbors
distable = table(round(patchcentersYX(:,2)),round(patchcentersYX(:,1)),NaN(npatches,1),NaN(npatches,1),'VariableNames',{'x','y','closest_neighbor','distpix'});
[idx dist] = knnsearch(patchcentersYX,patchcentersYX,'K',2);
distable.closest_neighbor = idx(:,2);
distable.distpix = dist(:,2);
distable.distum = [1/pixPerUm(zoom,scope)] * distable.distpix;
analysis.meandistum = mean(distable.distum);

% density recovery profile
patch_centers_image = sparse(false(size(patchimage)));
patchcentersind = sub2ind(size(patchimage),round(patchcentersYX(:,1)),round(patchcentersYX(:,2)));
patch_centers_image(patchcentersind) = true;
nintervals = (length(densrecov_intervals_um));
circtemplate = cell(nintervals,1);
densrecov_intervals_pix = round(pixPerUm(zoom,scope) * densrecov_intervals_um);
for i = 1:nintervals
    radpix = densrecov_intervals_pix(i);
    [xmesh ymesh] = meshgrid(0:2*radpix,0:2*radpix);
    circvals = sqrt([xmesh-radpix].^2 + [ymesh-radpix].^2);
    circtemplate{i} = circvals <= radpix;
end

densrecov_patchesPerSqmm = NaN(npatches,nintervals);
for indpatch = 1:npatches
    [y x] = deal(round(patchcentersYX(indpatch,1)), round(patchcentersYX(indpatch,2)));
    for indrad = 1:nintervals
        radpix = densrecov_intervals_pix(indrad);
        fullsize_plus_circle = false(size(patchimage));
        extenttop = max([1 y-radpix]);
        extentbot = min([size(patchimage,1) y+radpix]);
        extentleft = max([1 x-radpix]);
        extentright = min([size(patchimage,2) x+radpix]);
        templatecropped = circtemplate{indrad}(radpix-[y-extenttop]+1 : radpix+[extentbot-y]+1, radpix-[x-extentleft]+1 : radpix+[extentright-x]+1);
        fullsize_plus_circle(extenttop:extentbot,extentleft:extentright) = templatecropped;
        regiontosearch = roi & fullsize_plus_circle;
        searched_area_pix = numel(find(regiontosearch));
        searched_area_sqmm = searched_area_pix / [pixPerUm(zoom,scope)]^2 * 0.001^2;
        npatchesfound = length(find(patch_centers_image & regiontosearch))-1;
        densrecov_patchesPerSqmm(indpatch,indrad) = npatchesfound / searched_area_sqmm;
    end
end
distable.densrecov_patchesPerSqmm = densrecov_patchesPerSqmm;
    
npixroi = numel(find(roi)); 
analysis.densrecov_intervals_um = densrecov_intervals_um; 
analysis.roi_area_sqmm = npixroi/[pixPerUm(zoom,scope)]^2 * 0.001^2;
analysis.patches_per_sqmm = npatches / analysis.roi_area_sqmm;
analysis.patchimage = patchimage;
analysis.roi = roi;
analysis.patchesimage_roi = patchesimage_roi;
analysis.patch_centers_image = patch_centers_image;
analysis.npatches = npatches;
analysis.meandistpix = mean(distable.distpix);
analysis.zoom = zoom;
analysis.scope = scope;
analysis.distable = distable;


if plotting
    plot(analysis.densrecov_intervals_um,mean(densrecov_patchesPerSqmm))
    xlabel('Distance from patch center (um)')
    ylabel('Density of M2 peaks per sq. mm')
end