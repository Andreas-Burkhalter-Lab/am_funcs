function Bo = sngdecimate(B,navg)
% SNGDECIMATE: average nearby time bins of sonogram
% Syntax:
%   Bo = sngdecimate(B,navg)
% where
%   B is the abs(sonogram) (or it could be log(abs(B/const)))
%   navg is the number of successive columns to average
% and
%   Bo is the output matrix.
%
% size(B,2) must be an integer multiple of navg.
  
% Copyright 2001 by Timothy E. Holy <holy@pcg.wustl.edu>

  
  [nfreqs,ntimes] = size(B);
  nto = ntimes/navg;
  if (ceil(nto) > nto)
    error('size(B,2) is not an integer multiple of navg');
  end
  
  Bo = reshape(B,[nfreqs, navg, nto]);
  Bo = squeeze(sum(Bo,2));
