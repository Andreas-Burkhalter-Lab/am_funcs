function intens = imroi_holylab(ip,roi,options)
% IMROI_HOLYLAB: calculate intensities in ROIs
% Syntax:
%   intensity = imroi(ip,roi)
%   intensity = imroi(ip,roi,options)
% where
%   ip is an imphys structure array;
%   roi is a region description; if numeric, it assumes a 3-by-nregions
%     matrix, where the first two rows are x and y, respectively, and the
%     3rd the radius (all in pixels);
%   options is a structure array with the following fields:
%     showroi: if true, opens a figure and draws a circle at the position
%       of each ROI;
%     showimg: if present, uses this as the image for the ROI figure,
%       rather than the default (first frame);
%     showindx: if true, draws blobs rather than circles on the frame;
% and
%   intensity is a nframes-by-nregions matrix, giving the mean intensity in
%     each ROI as a function of frame.
  
% Copyright 2005 by Timothy E. Holy

% Name changed 2011 (TEH) to imroi_holylab to avoid name conflict with function
% in the image processing toolbox.
  
  nframes = length(ip);
  nregions = size(roi,2);
  
  if (nargin < 3)
    options = struct;
  end
  if ~isfield(options,'showroi')
    options.showroi = 0;
  end
  if ~isfield(options,'showindx')
    options.showindx = 0;
  end
  
  % Check to see if crop regions are all identical
  xrange = unique(reshape([ip.xrange],2,nframes)','rows');
  yrange = unique(reshape([ip.yrange],2,nframes)','rows');
  if (size(xrange,1) > 1 || ...
      size(yrange,1) > 1)
    error('Crop regions must be same in all frames');
  end
  
  % Calculate indices for ROIs
   for i = 1:nregions
      % Form the square region, then reject points
      r = roi(3,i);
      xr = round(roi(1,i)) + (-ceil(r):ceil(r));
      yr = round(roi(2,i)) + (-ceil(r):ceil(r));
      xlist = repmat(xr,1,length(xr));
      ylist = repmat(yr,length(yr),1);
      xylist = [xlist(:) ylist(:)]; % Here's the set of points within square of width r
      dxylist = xylist - repmat(roi(1:2,i)',size(xylist,1),1);
      distxylist = sqrt(sum(dxylist.^2,2));
      indxkeep = find(distxylist < r);
      % Convert x,y coords to an index
      sz = [diff(yrange)+1 diff(xrange)+1];
      imindx{i} = sub2ind(sz,xylist(indxkeep,2)-yrange(1)+1,...
                          xylist(indxkeep,1)-xrange(1)+1);
   end
   
   % If desired, show the figure
   if options.showroi
     if isfield(options,'showimg')
       img = options.showimg;
     else
       img = imphysfetch(ip(1));
     end
     if options.showindx
       for i = 1:nregions
         img(imindx{i}) = 0;
       end
       imagesc(img);
     else
       imagesc(img);
       colormap(gray)
       plotroi(roi);
     end
   end
         
   % Compute mean intensity inside ROIs in each frame
   intens = zeros(nframes,nregions);
   for i = 1:nframes
      img = imphysfetch(ip(i));
      fprintf('.');
      for j = 1:nregions
         intens(i,j) = mean(img(imindx{j}));
      end
   end
   fprintf('\n');
