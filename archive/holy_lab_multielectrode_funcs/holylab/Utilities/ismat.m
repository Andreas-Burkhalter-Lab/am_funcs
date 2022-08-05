function tf = ismat(filename)
% ismat: test whether a file is a matlab file
%
% Syntax:
%   tf = ismat(filename)
% where filename is the name of the file to analyze.

% Copyright 2012 by Timothy E. Holy

  [fid,msg] = fopen(filename);
  if (fid < 0)
    error(msg)
  end
  b = fread(fid,[6 1],'uint8');
  fclose(fid);
  tf = all(b == uint8('MATLAB')');
end
