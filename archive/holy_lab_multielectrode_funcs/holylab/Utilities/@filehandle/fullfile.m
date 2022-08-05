function fn = fullfile(fh)
% FILEHANDLE/FULLFILE: return full filename for a file handle
% Syntax:
%   filename = fullfile(fh)
  
% Copyright 2007 by Timothy E. Holy
  
  fn = fullfile(fh.abspathstr,fh.filename);
