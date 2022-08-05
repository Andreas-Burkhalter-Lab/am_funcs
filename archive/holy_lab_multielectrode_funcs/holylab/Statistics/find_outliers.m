function [outlierindx,d] = find_outliers(X,idx,c,options)
% FIND_OUTLIERS: find points too far from landmarks
% This is often useful in cleaning up problems from kmeans.
%
% Syntax:
%   outlierindx = find_outliers(X,idx,c)
%   outlierindx = find_outliers(X,idx,c,options)
%   [outlierindx,d] = find_outliers(...)
% where
%   X is your d-by-npts matrix of data points;
%   idx is the index of the landmark closest to each point (1-by-npts)
%   c is a d-by-nclusts matrix giving the landmark positions
% and
%   outlierindx is an index vector pointing to those X which lie too far
%     from the landmark;
%   d is the distance of each outlier from its associated landmark.
%
% The options may have the following fields:
%   threshfac (default 10): affects the definition of "too far" from each
%     landmark, the threshold is set at threshfac * median(distance).
%   plot (default 0): if true, plots each outlier waveform (red) and its
%     associated cluster center (black), and pauses for a keypress between
%     each pair. (Works only in 2d)
%
% See also: KMEANS_LAB.
  
% Copyright 2005 by Timothy E. Holy
  
  if (nargin < 4)
    options = struct;
  end
  if ~isfield(options,'threshfac')
    options.threshfac = 10;
  end
  if ~isfield(options,'plot')
    options.plot = 0;
  end

  clabel = agglabel(idx);
  nlabels = length(clabel);
  threshfac2 = options.threshfac^2; % mindist returns square-distance
  for i = 1:nlabels
    curpoints = clabel{i};
    if ~isempty(curpoints)
      dtmp = mindist(X(:,curpoints),c(:,i));
    end
    thresh = threshfac2*median(dtmp);
    outlierindx{i} = find(dtmp > thresh);
    d{i} = dtmp(outlierindx{i});
  end
  outlierindx = cat(2,outlierindx{:});
  d = sqrt(cat(2,d{:}));
  [outlierindx,p] = sort(outlierindx);
  d = d(p);
  if options.plot
    for i = 1:length(outlierindx)
      tindx = outlierindx(i);
      plot(X(:,tindx),'r')
      hold on
      plot(c(:,idx(tindx)),'k')
      hold off
      pause
    end
  end