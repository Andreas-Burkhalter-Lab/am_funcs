function [snip,time,header] = LoadSnip(filename,channel,options)
% LoadSnip: load the snippets from a given channel in one file
%   [snip,time,header] = LoadSnip(filename,channel,options)
% The third argument (options) is optional; a structure with the
% following fields:
%   timesonly: if present & true, the output is of the form
%     [time,header] rather than [snip,time,header]
%   machfmt:  allow files of different endian status to be opened;
%   maxsnip: at most maxsnip snippets will be read, with the default
%     behavior to read all the snippets
% snip is returned in units of volts. This is the only reliable choice for
% sorting, as the user might change the A/D scaling between files.
%
% See also: EPHYS.

if (nargin > 2)
  [header,fid] = readheader(filename,options);
else
  [header,fid] = readheader(filename);
end
width = header.sniprange(2)-header.sniprange(1)+1;
chindx = find(header.channels == channel);
if (isempty(chindx))
  snip = [];
  time = [];
  fclose(fid);
  return;
end
nsnips = header.numofsnips(chindx);
if (nargin > 2 & isfield(options,'maxsnip'))
  maxsnip = options.maxsnip;
  if (maxsnip >= 0 & maxsnip < nsnips)
    nsnips = maxsnip;
  end
end
fseek(fid,header.timesfpos(chindx),'bof');
time = fread(fid,nsnips,'int32');
header.thresh = header.thresh*header.scalemult + header.scaleoff;
if (nargin < 3 | ~isfield(options,'timesonly') | ~(options.timesonly))
  fseek(fid,header.snipsfpos(chindx),'bof');
  snip = fread(fid,[width,nsnips],'int16')*header.scalemult + ...
         header.scaleoff;
else
  snip = time;
  time = header;
end
fclose(fid);
