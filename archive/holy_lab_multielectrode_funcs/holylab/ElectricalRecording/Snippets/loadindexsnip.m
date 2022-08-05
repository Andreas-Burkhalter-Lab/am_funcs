function [snips,times] = loadindexsnip(filename,channels,index,options)
% LOADINDEXSNIP: read only a subset of snippets
% [snips,times] = loadindexsnip(filename,channels,index,options)
% where
%   filename is the name of the snippet file.
%   channels is a vector of channel numbers, or the string 'all'.
%   index gives the snippet numbers to be read. If there is more than one
%     channel being read, index must be a cell array of indices, one for
%     each channel. If only one channel is being read, then index may be
%     a vector.
%   options (optional) controls how the file is opened and read. If
%     options is a string, then it is taken to be the read precision
%     (default 'int16'). If it is a structure, then the two fields
%     "precision" and "machfmt" are tested. See help for fopen and
%     fread.
%     The field "tovolts," if present and true, will cause the snippets to
%     be converted to units of volts.
%
% If index is a vector (and only one channel is being read), then the output
% snips is a sniplen-by-nsnips matrix. If index is a cell array (several
% channels are being read), then snips is a cell array of such matrices.
% Similarly, times is a vector, or a cell array of vectors, depending on
% whether more than one channel is being read.
%
% See also: READINDEXSNIP, FOPEN, FREAD.

% Copyright 2001-2005 by Timothy E. Holy

  prec = 'int16';
  mach = 'n';
  tovolts = 0;
  if (nargin > 3)
    if ischar(options)
      prec = options;
    elseif isstruct(options)
      if isfield(options,'precision')
        prec = options.precision;
      end
      if isfield(options,'machfmt')
        mach = options.machfmt;
      end
      if isfield(options,'tovolts')
        tovolts = options.tovolts;
      end
    end
  else
    options = struct;
  end
  [h,fid] = readheader(filename,options);
  if ischar(channels)
    if strcmp(channels,'all')
      channels = h.channels;
    else
      error('Unrecognized channel string');
    end
  end
  [comchan,chani1,chani2] = intersect(channels,h.channels);
  if (length(comchan) < length(channels))
    error('Not all desired channels were recorded');
  end
  chanIndex = chani2(chani1);
  nchans = length(chanIndex);
  if (nchans > 1 | iscell(index))
    if (~iscell(index) | length(index) ~= nchans)
      error('Mismatch between index and number of channels');
    end
    snips = cell(1,nchans);
    for i = 1:nchans
      fseek(fid,h.snipsfpos(chanIndex(i)),'bof');
      snips{i} = readindexsnip(fid,index{i},h.sniprange,prec);
      if tovolts
        snips{i} = h.scalemult*snips{i} + h.scaleoff;
      end
      % If desired, read the times and subset them
      if (nargout > 1)
        fseek(fid,h.timesfpos(chanIndex(i)),'bof');
        ttemp = fread(fid,h.numofsnips(chanIndex(i)),'int32');
        times{i} = ttemp(index{i});
      end
    end
  else
    fseek(fid,h.snipsfpos(chanIndex),'bof');
    snips = readindexsnip(fid,index,h.sniprange,prec);
    if tovolts
      snips = h.scalemult*snips + h.scaleoff;
    end
    if (nargout > 1)
      fseek(fid,h.timesfpos(chanIndex),'bof');
      ttemp = fread(fid,h.numofsnips(chanIndex),'int32');
      times = ttemp(index);
    end
 end
end

  
