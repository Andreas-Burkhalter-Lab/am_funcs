function setspimg(h,varargin)
% SETSPIMG: sets properties (axis limits) for sparse image plots
% setspimg(h,field1,val1,...)
% where
%   h is the axis handle
%   field is the field name ('XLim' or 'YLim')
%   val is the value to use in setting the field
% This assumes that the sparse matrix and limit data are
% stored as a structure in the axis' 'UserData' field, using the
% following fields:
%    m holds the sparse matrix itself
%    xlim holds the 2-vector of x limits
%    ylim holds the 2-vector of y limits
% Other fields are optional:
%    climrank: if present, will set the upper colorscale limit to
%      saturation at the given rank. For example, setting this to 0.9
%      means that the brightest 10% of non-zero pixels will be saturated.
%    npix: if present, overrides the on-screen axis dimensions for the
%      purposes of creating a decimated image from the sparse image.
%    
% This function will rescale images via imagesc, so you should set
% the axis property 'NextPlot' to 'replacechildren'.
  spdata = get(h,'UserData');
  axu = get(h,'Units');
  set(h,'Units','pixels');
  axpos = get(h,'Position');
  set(h,'Units',axu);
  xlim = get(h,'XLim');
  ylim = get(h,'YLim');
  clim = [xlim;ylim];
  lim = clim;
  uselim = clim;
  tlim = [spdata.xlim;spdata.ylim];
  npix = axpos(3:4);
  if isfield(spdata,'npix')
    npix = spdata.npix;
  end
  outnpix = npix;
  % Compute the range of the matrix that will be used in making the
  % plot, and the target number of pixels occupied by that range.
  for i = 1:2:length(varargin)
    switch varargin{i}
     case 'XLim',
      index = 1;
     case 'YLim',
      index = 2;
     otherwise,
      error('Field name not recognized');
    end
    val = varargin{i+1};
    lim(index,:) = val;
    uselim(index,:) = IntersectIntervals(val,tlim(index,:));
    outnpix(index) = round(npix(index) * diff(uselim(index,:))/ ...
                           diff(val));
  end
  % Convert the range to index numbers and get the sub-matrix
  [ny,nx] = size(spdata.m);
  scalefac = diag(([nx ny] - [1 1])./diff(tlim'),0);
  subrange = round(scalefac * (uselim - repmat(tlim(:,1),1,2))) + 1;
  msub = spdata.m(subrange(2,1):subrange(2,2), ...
           subrange(1,1):subrange(1,2));
  % Decimate the matrix
  if any(size(msub) > outnpix)
    [i,j,s] = find(msub);
    stretchfac = (outnpix(2:-1:1) - [1 1]) ./ size(msub);
    if (stretchfac(1) < 1)
      i = round((i-1) * stretchfac(1) + 1);
    end
    if (stretchfac(2) < 1)
      j = round((j-1) * stretchfac(2) + 1);
    end
    msub = sparse(i,j,s,min(outnpix(2),size(msub,1)), ...
                        min(outnpix(1),size(msub,2)));
    % Normalize, so intensity is constant independent of binning
    % This was commented out because for a sparse matrix, the best way to
    % do this might well be not to normalize!
    %stretchfac = min([1 1;stretchfac]);
    %msub = msub*prod(stretchfac);
    %inz = find(msub);
    %hcfig = gcf;
    %figure
    %hist(msub(inz),30)
    %set(gca,'YScale','log')
    %figure(hcfig);
  end
  % Render the image
  imagesc(uselim(1,:),uselim(2,:),msub,'Parent',h);
  set(h,'XLim',lim(1,:),'YLim',lim(2,:));
  if isfield(spdata,'climrank')
    inz = find(msub);
    climrank = max(0,min(spdata.climrank,1));
    if (climrank ~= spdata.climrank)
      warning('climrank clipped to [0,1]');
    end
    indx = round(climrank*length(inz));
    if (indx > 0)
      snz = sort(msub(inz));
      set(h,'CLim',[0 snz(indx)]);
    end
  end
