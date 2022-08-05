function haxo = imsubplot(nrows,ncols,imsize,gap)
% IMSUBPLOT: create axes for image subplots
% This function insures that the figure and axis geometry use screen
% space efficiently.
%
% Syntax:
%   hax = imsubplot(nrows,ncols,imsize)
% where
%   nrows is the number of rows of images;
%   ncols is the number of columns of images;
%   imsize is the size of an image, imsize = size(im) (identical for each
%     image)
% 
% If you want to leave a larger gap at the top of each image (e.g., for a
% title), use the syntax
%   hax = imsubplot(nrows,ncols,imsize,gap)
% where gap (default 0.05) is the fraction of the image height used to
% separate the images.
  
% Copyright 2005 by Timothy E. Holy
  
  clf   % Clear the current figure
  set(gcf,'MenuBar','none')
  set(gca,'Position',[0.02 0.02 0.96 0.96]); % Use most of the space
  if (nargin < 4)
    gap = [];
  end

  % Arrange the subplots
  splits = isvdf_axsplits(nrows,gap);
  if ~isempty(splits)
    haxv = SplitVert(splits);
    if isempty(gap)
      delete(haxv(2:2:end));
      haxv = haxv(1:2:end);
    else
      delete(haxv(1:2:end));
      haxv = haxv(2:2:end);
    end
  else
    haxv = gca;
    if ~isempty(gap)
      set(gca,'Position',[0.02 0.02 0.96 0.96-gap]); % Use most of the space
    end
  end
  for j = 1:nrows
    if (ncols > 1)
      splits = isvdf_axsplits(ncols);
      haxtmp = SplitHoriz(splits,haxv(j));
      delete(haxtmp(2:2:end));
      haxtmp = haxtmp(1:2:end);
      hax(:,j) = haxtmp(:);
    else
      hax(:,j) = haxv(j);
    end
  end

  % Set up axes for image plots
  axis(hax(:),'image');
  set(hax(:),'Visible','off')

  % Set the figure size appropriately
  npix = imsize.*[nrows ncols];
  set(gcf,'Units','normalized');
  set(gcf,'Position',[0 0 npix([2 1])/max(npix)/2]);

  % If desired, output the axis handles
  if (nargout > 0)
    haxo = hax;
  end
  
function splits = isvdf_axsplits(n,gap)
  if (nargin < 2 || isempty(gap))
    gap = 0.05;
    last = 2*n-1;
  else
    last = 2*n;
  end
  splits = repmat([1 gap],1,n);
  splits = cumsum(splits(1:last));
  splits = splits/splits(end);
  if length(splits > 1)
    splits = splits(1:end-1);
  else
    splits = [];
  end
