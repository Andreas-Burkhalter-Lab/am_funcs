function masks=mask_objs(subframe, nMice, isShowInterResult, textureFilter, ch, offset)
% todo: maybe I should pick another func name
% todo: change nMice to nObj

% entrofilt first then convert to gray so we can take color info
%    which is better than "convert to gray then do entrofilt()"

% smooth then entropyfilt
binsize=2; % 25  % 10
% image_enfilt = entropyfilt(subframe/binsize); % NOTE: subframe is uint8, so /binsize is kind of smoothing
% image_enfilt = rangefilt(subframe/binsize); % NOTE: subframe is uint8, so /binsize is kind of smoothing
% stdfilt() is way fast but the result is too bad
% image_enfilt = stdfilt(subframe/binsize); % NOTE: subframe is uint8, so /binsize is kind of smoothing
chIndices={1, 2, 3, 1:3};
idx=strmatch(ch, {'r', 'g', 'b', 'a'}); % 'exact' is not used for things like 'all', 'red' etc.
if(isempty(idx))
   error('incorrect ch param');
else
   ch=chIndices{idx};
end

image_enfilt = textureFilter(subframe(:,:,ch)/binsize); % NOTE: subframe is uint8, so /binsize is kind of smoothing
if isShowInterResult
   fold=floor(255/max(image_enfilt(:)));
   figure; image(uint8(fold*image_enfilt)); axis image; title('image_enfilt', 'interpr', 'none')
end

% convert to gray
if(idx==4)
   image_enfilt_gray=color2gray(image_enfilt);
else
   image_enfilt_gray=image_enfilt;
end
if isShowInterResult
   figure; colormap(gray(256)); imagesc(image_enfilt_gray);
   axis image; title('image_enfilt_gray', 'interpr', 'none')
end

% convert to 8 bit and invert
image_enfilt_gray_8_inv=255-to8bit(image_enfilt_gray);
if isShowInterResult
   figure; colormap(gray(256));
   image(image_enfilt_gray_8_inv) % after make it in range 0-255, then no need call imagesc()
   axis image; title('image_enfilt_gray_8_inv', 'interpr', 'none')
end


% thresholding and convert to black-white image (0 or 255)
[n,x]=hist(image_enfilt_gray_8_inv(:), 64);
[y,i]=max(n);
peak=round(x(i));
thresh=(255-peak)*offset+peak;

image_enfilt_gray_8_inv_thresh=image_enfilt_gray_8_inv;
image_enfilt_gray_8_inv_thresh(image_enfilt_gray_8_inv_thresh<thresh)=0;
image_enfilt_gray_8_inv_thresh(image_enfilt_gray_8_inv_thresh~=0)=255;
if isShowInterResult
   figure; colormap(gray(256)); image(image_enfilt_gray_8_inv_thresh),
   axis image; title('image_enfilt_gray_8_inv_thresh', 'interpr', 'none')
end

bw=image_enfilt_gray_8_inv_thresh;

% % clear border
% image_enfilt_gray_8_inv_thresh_clearborder=imclearborder(image_enfilt_gray_8_inv_thresh, 4);
% if isShowInterResult
%    figure; colormap(gray(256));
%    imagesc(double(image_enfilt_gray_8_inv_thresh_clearborder)),
%    axis image; title('image_enfilt_gray_8_inv_thresh_clearborder', 'interpr', 'none')
% end

% % now convert to black&white image (0 or 255)
% % 
% % bw=image_enfilt_gray_8_inv_thresh;
% bw=image_enfilt_gray_8_inv_thresh_clearborder;
% 
% bw(bw~=0)=255;
% if isShowInterResult
%    figure; colormap(gray(256)); image(double(bw)),
%    axis image; title('bw', 'interpr', 'none')
% end

% for stdfilt(), we need clean up the details a bit. Note that a lowpass
%   filter doens't help much though.
% Also medfilt2() doesn't help either
if(isequal(textureFilter, @stdfilt))
   se = strel('disk', 4);
   imopen_bw=imopen(bw, se);
   if isShowInterResult
      figure; colormap(gray(256)); image(double(imopen_bw));
      axis image; title('imopen_bw', 'interpr', 'none')
   end
   bw_old=bw; % for debug only
   bw=imopen_bw;
end

% gaussian filter to remove small dots then round to black-white again
bw_gau=imfilter_gaussian(single(bw), [2 2]);
bw_gau_bw=bw_gau;
bw_gau_bw(bw_gau_bw<128)=0;
bw_gau_bw(bw_gau_bw>=128)=255;
if isShowInterResult
   figure; colormap(gray(256));
   image(double(bw_gau_bw)), axis image; title('bw_gau_bw', 'interpr', 'none')
end

% imfill, bwselect, but bwlabel is the best
objLabel=bwlabel(bw_gau_bw);
if isShowInterResult
   figure; rgb=label2rgb(objLabel); imshow(rgb); axis image; title('labeled objs')
end

% now create the masks for the mice
nObj=max(objLabel(:));
masks=cell(1,nObj);
for objIndex=1:nObj
   masks{objIndex}=objLabel==objIndex;
end
% sort masks by area sizes
areas=cellfun(@(x)sum(x(:)), masks);
[areaSizes, indices]=sort(areas, 'descend');
masks=masks(indices);
textureMasks=masks;
% now clean up masks by area size and grow masks by intensity
ch=subframe(:,:,2); % TODO: the channel to use
intenMasks={};
for maskIndex=1:length(masks)
   mask=masks{maskIndex};
   if(isShowInterResult)
      figure; rgb=label2rgb(mask); imshow(rgb); axis image; title(['texture mask #' num2str(maskIndex)])
   end
   if(areaSizes(maskIndex)<81 && length(intenMasks)>=nMice) % NOTE: hard coded 81
      break;
   end
   masked=double(ch(mask));
   m=mean(masked); s=std(masked);
   if(s/m<0.01 && m>255*0.95) % NOTE: hard coded
      % an overexposed area
      continue;
   end
   if(s/m>1) % NOTE: hard coded
      % don't know what to do, maybe an area covering two mice
      intenMasks{end+1}=mask;
      continue;
   end
   intenMasks{end+1}=maskFromInten(mask, ch, isShowInterResult);
   if(isShowInterResult)
      figure; rgb=label2rgb(intenMasks{end}); imshow(rgb); axis image; title(['inten mask #' num2str(length(intenMasks))])
   end
end % for, each mask based on texture

masks=intenMasks;
% masks=textureMasks;

% get rid of masks that are way too narrow
indices=[];
for idx=1:length(masks)
   mask=masks{idx};
   [rows, cols]=find(mask);
   diffRow=max(rows)-min(rows);
   diffCol=max(cols)-min(cols);
   ratio=diffRow/diffCol;
   if(ratio<1)
      ratio=1/ratio;
   end
   if(ratio>5) % TODO: hard coded
      indices(end+1)=idx;
   end
end
masks(indices)=[];

if(length(masks)<nMice)
   return; % don't know what to do, and reuse previous result
end

% sort masks again

areas=cellfun(@(x)sum(x(:)), masks);
[areaSizes, indices]=sort(areas, 'descend');
masks=masks(indices(1:nMice));

% now we have to deal w/ cases when mice's masks are connected, that is, to
% split connected masks. Here we use the intensities to separate them.
image_gray=color2gray(subframe);
countPerMask=round(areaSizes(1:nMice)/mean(areaSizes(1:nMice)));
if(isequal(countPerMask, [1 1 0]))
   countPerMask=[1 1 1];
end
masksOut={};
for idx=1:nMice
   if(countPerMask(idx)==0) continue; end
   mask=masks{idx};
   if(countPerMask(idx)==1) masksOut{end+1}=mask; continue; end
   % if a mask covers two mice, actually using the mean (b/c the areas
   % roughtly equal) to split is much fast. But here I use kmean for more general
   % processing.
   masked=image_gray(mask); % here masked area is vectorized
   tLabels=kmeans(masked, countPerMask(idx));
   % convert the labeling back to matrix
   indices=find(mask);
   labels=double(mask);
   labels(indices)=tLabels;
   % ok, now the splitted masks
   splited={};
   for newMaskIdx=1:countPerMask(idx)
      splited{end+1}=labels==newMaskIdx;
      % masksOut{end+1}=labels==newMaskIdx;
   end
   splitedAreas=cellfun(@(x)sum(x(:)), splited);
   [splitedAreaSizes, tt]=sort(splitedAreas, 'descend');
   splitedCom=cellfun(@centerOfMass, splited, 'UniformOutput', false);
   distance=sum((splitedCom{1}-splitedCom{2}).^2)^0.5;
   if(splitedAreaSizes(1)/splitedAreaSizes(2)>7 ) || (distance<10) % TODO: hardcoded 10 pixels (about 2cm)
      % a failed splition
      masksOut{end+1}=mask;
      zeroIndices=find(countPerMask==0);
      countPerMask(zeroIndices(1:min(countPerMask(idx)-1, length(zeroIndices))))=1;
   else
      masksOut=[masksOut splited];
   end
end % for, each original mask
masks=masksOut;

if isShowInterResult
   tLabel=masks{1}; % for plotting only
   for mouseIndex=2:nMice
      tLabel=tLabel+masks{mouseIndex}*mouseIndex;
   end
   figure; tRgb=label2rgb(tLabel); imshow(tRgb); axis image; title('found mice')
end


