function header = resolvelegacyheader(fid,magicnum,options)
% RESOLVELEGACYHEADER: a gateway which deals with old header types
%
% See also: READHEADER
  
  if (max(abs(magicnum(1:4) - [196   228    43    59])) == 0)
    % This is a WashU v1 style header
    fseek(fid,0,'bof');
    header = ReadAIHeaderWU1(fid); 
    % Note: could check the byte-swapped version of the above magic
    % number to yield an automatic method; not yet implemented
  else
    % Assume this is an ancient Harvard-style header
    % Determine whether it's a snippet file, an envelope file, etc.
    fseek(fid,0,'bof');
    headersize = fread(fid,1,'uint32');
    type = fread(fid,1,'int16');
    fseek(fid,0,'bof');
    if (type == 1)
      if (strcmp(options.headertype,'ai'))
        header = ReadAIHeaderHarvard(fid);
      elseif (strcmp(options.headertype,'snip'))
        header = ReadSnipHeader(fid);
      else
        error('Need to know header type, ai or snip');
      end
    elseif (type == 4)
      header = ReadEnvHeader(fid);
    elseif (type == 8)
      header = ReadVidHeader(fid);  % What function should this be,
                                    % anyway?
    else
      error('Unrecognized legacy header type');
    end
  end
