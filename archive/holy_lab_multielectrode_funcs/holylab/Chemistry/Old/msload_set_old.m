function [Mc,mz,tc] = msload_set(filenames,mzoffset,tShift,tRange,options)
% MSLOAD_SET: load a set of mass spec experiments to matrix format
% Syntax:
%   [Mc,mz,tc] = msload_set(filenames)
%   [Mc,mz,tc] = msload_set(filenames,mzoffset)
%   [Mc,mz,tc] = msload_set(filenames,mzoffset,tShift,tRange)
%   [Mc,mz,tc] = msload_set(filenames,mzoffset,tShift,tRange,options)
% where
%   filenames is a cell array of mzXML file names;
%   mzoffset is a 1-by-n_files vector specifying the offset to apply to the
%     m/z index when loading the matrices (this is used to perform
%     "m/z registration")
%   tShift is a 1-by-n_files vector specifying the temporal shift (in
%     units of minutes) to apply to each set of scan times
%   trange is a 1-by-2 vector listing the range of scan times to load and
%     convert to matrix form
%   options is a structure that can be used to control the behavior of
%     scan2matrix and mz2indx
% and
%   Mc is a cell array of sparse matrices, as loaded by scan2matrix;
%   mz is a vector of m/z values corresponding to each row of Mc{i};
%   tc is a cell array of times (in minutes) associated with each scan in
%     each file.
%
% See also: SCAN2MATRIX.

% Copyright 2009 by Timothy E. Holy

  if ~iscell(filenames)
    filenames = {filenames};
  end
  n_files = length(filenames);
  if (nargin < 2)
    mzoffset = zeros(1,n_files);
  end
  if any(mzoffset < 0)
    error('All mzoffset must be nonnegative; subtract off the smallest value');
  end
  if (nargin < 5)
    options = struct;
  end

  tc = cell(1,n_files);
  Mc = cell(1,n_files);
  for fileIndex = 1:n_files
    % Load the data from a single file and keep just the profile data
    fprintf('Loading %s...\n',filenames{fileIndex});
    scan = msload_mzXML(filenames{fileIndex});
    flag = [scan.ms_level] == 1;
    scan = scan(flag);
    st = [scan.scan_time]/60;
    tc{fileIndex} = st;
    if (nargin > 2)
      st = st + tShift(fileIndex);
      flagst = st > tRange(1) & st < tRange(2);
      scan = scan(flagst);
      tc{fileIndex} = st(flagst);
    end
    if (fileIndex == 1)
      % If it's the first file, calculate the conversion between m/z and row
      % numbers
      [Mc{fileIndex},f,finv] = scan2matrix(scan,[],mzoffset(fileIndex),options);
      mxmzI = ceil(f(max([scan.mz]))+max(mzoffset));
    else
      % Use the same conversion used for the first file
      Mc{fileIndex} = scan2matrix(scan,f,mzoffset(fileIndex));
    end
  end
  % Make sure they all have the same size
  for i = 1:length(Mc)
    if (size(Mc{i},1) < mxmzI)
      Mc{i}(mxmzI,end) = 0;
    end
  end
%   M = cat(2,Mc{:});
%   t = cat(2,tc{:});
  mz = finv(1:mxmzI)';
  fprintf('Done\n');
end
