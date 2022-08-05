function [M,mz2i,i2mz] = scan2matrix(scan,mz2i,offset,options)
% SCAN2MATRIX: convert mass spec data to sparse matrix format
% Syntax:
%  M = scan2matrix(scan)
% This is the most basic syntax. Here, scan is a structure of individual
% scans, such as returned by msload_mxXML. M is the matrix, with the index
% corresponding to m/z as the first index and time as the second.
%
%  [M,mz2i,i2mz] = scan2matrix(scan)
% This returns the functions for converting between m/z and the
% corresponding index. See mz2indx. You can get the m/z value corresponding
% to each index with
%    mz = i2mz(1:size(M,1));
%
%  M = scan2matrix(scan,mz2i)
% This uses a predefined function for converting m/z to an index. This is
% useful if you are converting multiple runs and want to use a consistent
% coordinate.
%
% [M,mz2i,i2mz] = scan2matrix(scan,[],offset)
% [M,...] = scan2matrix(scan,mz2i,offset)
% Use these forms if you have performed "m/z registration" to align multiple
% runs to compensate for drifting m/z values. Note the offset should be
% non-negative (so subtract off the minimum offset across files).
%
% [M,mz2i,i2mz] = scan2matrix(scan,[],0,options)
% Use this form if you want to pass options to mz2int.
%
% See also: MZ2INDX.


% Copyright 2009 by Timothy E. Holy

  if (nargin < 4)
    options = struct;
  end
  n_scans = length(scan);
  min_mz = min([scan.mz]);
  if (nargin < 2 || isempty(mz2i))
    % Figure out how to convert m/z values to an index
    tI = [scan.totIntensity];
    [mxI,maxIndex] = max(tI);
    [mz2i,i2mz] = mz2indx(scan(maxIndex).mz,min_mz,options);
  end
  if (nargin < 3)
    offset = 0;
  end
  % Convert each scan
  mzI = cell(1,n_scans);
  sI = mzI;
  for i = 1:n_scans
    mzI{i} = mz2i(scan(i).mz);
    sI{i} = repmat(i,1,length(mzI{i}));
  end
  % Use linear interpolation to fill the intensities
  mzIndex = double(cat(2,mzI{:})) + offset;
  mzInt = floor(mzIndex);
  frac = mzIndex - mzInt;
  scanIndex = cat(2,sI{:});
  I = double([scan.intensity]);
  M = sparse([mzInt, mzInt+1],[scanIndex, scanIndex],[I.*(1-frac), I.*frac]);
end