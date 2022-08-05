function d = whisdist(w1,w2,options)
% WHISDIST: distance between chirps using warping
% Syntax:
%   d = whisdist(w1,w2)
%   d = whisdist(w1,w2,options)
% where
%   w1 and w2 are the pitches of two chirps;
%   options may have the following fields:
%     freeends: if false, aligns by DTW, whereas if true, uses
%       DTWFE (default false);
%     frac: if freeends is false, is computes the error using only the
%       frac best-aligned points (default 0.95);
% and
%   d is the distance between them.
%
% See also: DTW,DTWFE.
  
% Copyright 2005 by Timothy E. Holy
  
  if (nargin < 3)
    options = struct;
  end
  if ~isfield(options,'freeends')
    options.freeends = 0;
  end
  if ~isfield(options,'frac')
    options.frac = 0.95;
  end
  
  if options.freeends
    [dist,mapping,ww1,ww2,ellmin] = dtwfe(w1,w2,struct('ellpower',2,'minoverlap',3));
    d = mean(abs(ww1-ww2));
    % Adjust for shifts
    corrfac = (length(w1) + length(w2)) / (2*ellmin);
    d = d * corrfac;
  else
    [dist,mapping,ww1,ww2] = dtw(w1,w2);
    d = sort(abs(ww1-ww2));
    d = mean(d(1:ceil(options.frac*length(ww1))));
  end
