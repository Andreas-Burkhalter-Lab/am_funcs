function [Mc,mz2i,i2mz,tc] = msprof_load_set(filenames,options)
% MSPROF_LOAD_SET: load a set of mass spec experiments to matrix format
% Syntax:
%   [Mc,mz2i,i2mz,tc] = msprof_load_set(filenames)
%   [Mc,mz2i,i2mz,tc] = msprof_load_set(filenames,options)
% where
%   filenames is a cell array of mzXML file names;
%   options is a structure that may have the following fields:
%     mzi_shift is a 1-by-n_files vector specifying the offset to apply to the
%       m/z index when loading the matrices (this is used to perform
%       "m/z registration")
%     t_shift is a 1-by-n_files vector specifying the temporal shift (in
%       units of minutes) to apply to each set of scan times
%     trange is a 1-by-2 vector listing the range of scan times to load and
%       convert to matrix form, in minutes (note this is different from
%       msload_mzXML, which is in seconds)
% and
%   Mc is a cell array of "sparse matrices" parametrized in structure-format;
%   mz2i is a function handle converting m/z values to indices;
%   i2mz is a function handles converting indices to m/z values;
%   tc is a cell array of scan times, in minutes (not seconds like mzXML).
%
% See also: SCAN2MATRIX.

% Copyright 2009 by Timothy E. Holy

  if ~iscell(filenames)
    filenames = {filenames};
  end
  n_files = length(filenames);
  if (nargin < 2)
    options = struct;
  end
  options = default(options,...
    'mzi_shift',zeros(1,n_files),...
    't_shift',zeros(1,n_files),...
    'trange',[0 Inf]);
  if any(options.mzi_shift < 0)
    error('All mzi_shift must be nonnegative; subtract off the smallest value');
  end

  tc = cell(1,n_files);
  Mc = cell(1,n_files);  
  for fileIndex = 1:n_files
    % Load the data from a single file and keep just the profile data
    fprintf('Loading %s...\n',filenames{fileIndex});
    scan = msload_mzXML(filenames{fileIndex},struct('ms_level',1, ...
      'scan_time_range',60*(options.trange - options.t_shift(fileIndex))));
    st = [scan.scan_time]/60 + options.t_shift(fileIndex);
    tc{fileIndex} = st;
    if (fileIndex == 1)
      % If it's the first file, calculate the conversion between m/z and row
      % numbers, using the scan with highest intensity as a study case
      min_mz = min([scan.mz]);
      tI = [scan.totIntensity];
      [mxI,maxIndex] = max(tI);
      [mz2i,i2mz] = mz2indx(double(scan(maxIndex).mz),double(min_mz),options);
    end
    n_scans = length(scan);
    % Initialize each "sparse matrix" in column-major form
    Mc{fileIndex} = repmat(struct('rowi',[],'value',[],'rowi_min',[],'rowi_max',[]),1,n_scans);
    for scanIndex = 1:n_scans
      I = scan(scanIndex).intensity;
      keepflag = I > 0;
      I = double(I(keepflag));
      mz = double(scan(scanIndex).mz(keepflag));
      % Convert to an index
      mzIndex = mz2i(mz) + options.mzi_shift(fileIndex);
      % Convert indices to integers, assigning fractional components to
      % the integer on either side. We go to some effort to avoid just
      % duplicating the array and then calling sort, because sort is slow.
      mzInt = floor(mzIndex);
      frac = mzIndex - mzInt;
      mzIntm = [mzInt; mzInt+1];
      vm = [(1-frac).*I; frac.* I];
      indexmerge = find(diff(mzInt) == 0);
      vm(:,indexmerge) = vm(:,indexmerge) + vm(:,indexmerge+1);
      mzIntm(:,indexmerge+1) = [];
      vm(:,indexmerge+1) = [];
      indexmerge = find(diff(mzIntm(1,:)) == 1);
      vm(1,indexmerge+1) = vm(1,indexmerge+1) + vm(2,indexmerge);
      vm(2,indexmerge) = 0;
      mzInt = mzIntm(:);
      v = vm(:);
      keepflag = v>0;
      mzInt = mzInt(keepflag);
      v = v(keepflag);
      Mc{fileIndex}(scanIndex).rowi = mzInt;
      Mc{fileIndex}(scanIndex).value = v;
      Mc{fileIndex}(scanIndex).rowi_min = mzInt(1);
      Mc{fileIndex}(scanIndex).rowi_max = mzInt(end);
    end
  end
  fprintf('Done\n');
end
