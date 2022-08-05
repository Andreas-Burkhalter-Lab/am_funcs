function plotroi(roidata,roioptions)
%PLOTROI: plot ROIs
%Syntax:
%  plotroi(roidata)
%  plotroi(roidata,roioptions)
%where
%  roidata specifies the ROI geometry; if it's numeric, it's assumed
%    these are circles specified by a 3-by-nrois matrix, where rows 1 and
%    2 give the x- and y-coordinates of the center of each circle, and
%    row 3 gives the radius;
%  roioptions is an optional structure with the following fields:
%    color (optional): a #ROIs-by-3 matrix specifying the color used for
%      each ROI;
%    linewidth (optional): a scalar (to use the same value for each ROI)
%      or vector (to specify a particular value for each ROI) setting the
%      linewidth;

% Copyright 2005 by Timothy E. Holy
  
  if (nargin < 2)
    roioptions = struct;
  end
  if ~isfield(roioptions,'color')
    roioptions.color = get(gca,'ColorOrder');
  end
  if ~isfield(roioptions,'linewidth')
    roioptions.linewidth = 2;
  end
  
  if isnumeric(roidata)
    nregions = size(roidata,2);
    if isscalar(roioptions.linewidth)
      roioptions.linewidth = repmat(roioptions.linewidth,1,nregions);
    end
    theta = linspace(0,2*pi,500);
    co = roioptions.color;
    for i = 1:nregions
      xcirc = roidata(1,i) + roidata(3,i)*cos(theta);
      ycirc = roidata(2,i) + roidata(3,i)*sin(theta);
      line(xcirc,ycirc,'Color',co(mod(i-1,size(co,1))+1,:),...
           'LineWidth',roioptions.linewidth(i));
    end
  else
    error('No other geometry implemented yet');
  end
