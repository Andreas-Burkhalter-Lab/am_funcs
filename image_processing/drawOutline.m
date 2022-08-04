% draws a border around all shapes within the input image of specified
% thickness
%[outlinePix,bounds,labelmat] = drawOutline(img,thickness_pix)
%%% outlinePix_yx = cell containing coordinate subcripts of border of specified thickness for each shape in input image
%%% innerPix_yx = cell containing coordinate subcripts of pixels within each shape, including the border 
%%% bounds_yx = boundaries of all shapes (inside and excluding outlinePix)
%%% labelmat = image of shapes in input image, with different values for each shape
% img = input black/white image; shapes = 1, background = 0
% thickness_pix = thickness of border to be drawn around each shape
% borderCannotOverlapShapes = if == 'noOverlap', eliminate from border regions all pixels contained by a shape  
%%% last updated 3/16/18

function [outlinePix_yx,innerPix_yx,bounds_yx,labelmat] = drawOutline(img,thickness_pix,borderShapeOverlap)

if ~exist('borderShapeOverlap','var')
    noBorderShapeOverlap = false;
elseif strcmp(borderShapeOverlap,'noOverlap')
    noBorderShapeOverlap = true; % eliminate from border regions all pixels contained by a shape  
end
[bounds_yx,labelmat] = bwboundaries(img,'noholes'); % 'noholes' to not treat FALSE pixels within patches as patch pixels
npatches = length(bounds_yx);
outlinePix_yx = cell(npatches,1);
innerPix_yx = cell(npatches,1);
maxY = size(img,1);
maxX = size(img,2);
shapePixelInds = find(img);

for patchind = 1:npatches
    clear intpatchSubsYX
    % get innerPix
    thisPatchImage = labelmat == patchind;
    [innerPix_yx{patchind}(:,1) innerPix_yx{patchind}(:,2)] = ind2sub(size(thisPatchImage),find(thisPatchImage)); % store as x,y from top,left
    
    % draw border
    [intpatchSubsYX(:,1) intpatchSubsYX(:,2)] = find(bwdist(thisPatchImage) < thickness_pix & bwdist(thisPatchImage) > 0);
    if noBorderShapeOverlap % eliminate from border regions all pixels contained by a shape  
        intpatchinds = sub2ind(size(img),round(intpatchSubsYX(:,1)),round(intpatchSubsYX(:,2)));
        [~,overlapPixels] = intersect(intpatchinds,shapePixelInds); %get index of pixels within intpatchinds that overlap shape pixels
        intpatchinds(overlapPixels) = []; % delete overlapping pixels
        [outlinePix_yx{patchind}(:,1),outlinePix_yx{patchind}(:,2)] = ind2sub(size(img),intpatchinds);
    else
        outlinePix_yx{patchind} = round(intpatchSubsYX);
    end
end