function [fh,msgo] = fopen(fh,permission)
% FILEHANDLE/FOPEN
% Syntax:
%   fh = fopen(fh)
% Opens for reading
%   fh = fopen(fh,permission)
% Opens in reading, writing, etc. modes. See help for the builtin FOPEN.
%   [fh,msg] = fopen(...)
% Returns a message about the success or failure.
%
% Note when you're reading info you'll still want to query uselfs to see
% which fread utilities to use.  
%
% See also: FOPEN.
  
  if (nargin < 2)
    permission = 'r';
  end
  
  % Test for incompatibilities in machineformat
  [c,maxsize,machineformat] = computer;
  if (fh.uselfs && isempty(strmatch(lower(fh.machfmt),...
                                    {'n',lower(machineformat)})))
    error('uselfs works only with native ordering');
  end
  
  fullpath = fullfile(fh.abspathstr,fh.filename);

  if fh.uselfs
    %[fh.fid,msg] = openlfs(fullpath,permission);
    [fh.fid,msg] = openlfs(fullpath);
  else
    [fh.fid,msg] = fopen(fullpath,permission,fh.machfmt);
  end
  
  if (nargout > 1)
    msgo = msg;
  end
  