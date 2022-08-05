function [mz,ic] = msprof_parse(filename)
% MSPROF_PARSE: parse profiled MS data
% Syntax:
%   [mz,ic] = msprof_parse(filename)
% where
%   filename is the name of the file
% and
%   mz is a vector containing the m/z values scanned;
%   ic is the ion current, stored as a matrix ic(m/z index,sweep #).
%     Sweeps occur over time as the sample progresses into the MS.
%
% See also: MSPROF_FACTOR.

  [fid,msg] = fopen(filename,'r');
  if (fid < 0)
    error(msg);
  end
  A = fscanf(fid,'%g');
  mz = A(2:2:end);  % m/z
  ic = A(1:2:end);  % ion current
  % Find the # of m/z values in each sweep
  mzlen = find(diff(mz) < 0,1,'first');
  mz = mz(1:mzlen);
  nscans = floor(numel(ic)/mzlen);
  if (nscans * mzlen < numel(ic))
    warning('Partial scan detected, truncating');
  end
  ic = reshape(ic(1:mzlen*nscans),mzlen,nscans);
  