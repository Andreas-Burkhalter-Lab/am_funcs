function sto = fclose(fh)
% FILEHANDLE/FCLOSE
% Syntax:
%   fclose(fh)
%   st = fclose(fh) also works, if you're not using Large File System
%     utilities (this could be fixed if this syntax is necessary)
%
% See also: FCLOSE.
  
  if fh.uselfs
    closelfs(fh.fid);
  else
    st = fclose(fh.fid);
    if (nargout > 0)
      sto = st;
    elseif (st < 0)
      % If the user isn't collecting the output, then issue an error if
      % there was a problem
      error(['Closing file ' fh.filename ' was not successful']);
    end
  end
  