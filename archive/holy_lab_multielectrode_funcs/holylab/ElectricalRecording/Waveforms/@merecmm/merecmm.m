function memm = merecmm(filename,varargin)
% MERECMM: memory map a MEREC file
% This gives convenient access to the parsed header as well as the data
% of a MEREC file.  The main advantage of this is simplified access for the
% programmer, without sacrificing reasonable performance. However, there is
% some overhead; you can go faster yet if you directly control the reading.
%
% Syntax:
%   memm = merecmm(filename)
%   memm = merecmm(filename,param1,value1,...)
% where filename is the name of the .merec file you want access to.
%
% Once you've done this, you can access data in the following way:
%   memm.contiguous = false;
%   x = memm([0 37 52],10000:50:1e7)
% would load data for channels 0, 37, and 52 (note these are channel #s,
% not channel indices!) for every 50th scan, from scan # 10000 to the
% 10-millionth.
% You can also do this:
%   memm.contiguous = true;
%   x = memm([0 37 52],[10000 15000]);
% would load data for channels 0, 37, and 52 (note these are channel #s,
% not channel indices!) for scans # from 10000 to 15000. The "contiguous"
% parameter determines whether the second index is interpreted as the
% limits of a range (when contiguous is true) or to be interpreted
% literally (when contiguous is false).
%
% A possible future way to access the data is with a single index:
%   x = memm(4:64:60000)
% would read all the data from the 4th channel (note here 4 is a channel
% index, not channel # as above!), if 64 channels were recorded. In terms
% of scans, it would include the first up to the 60000th A/D conversion
% (which occuring during the 938th scan), thus returning a vector of length
% 938.
%
% You can also get auxillary data like the following:
%   memm.header: returns the parsed MEREC header (like READHEADER)
%   memm.channels: returns the channel #s
%   memm.scanrate: returns the # of scans per second
%   memm.nscans: returns the total # of scans in the file
%   memm.size: returns the size of the data array
%   memm.filename: returns the name of the mapped file
%   memm.tovolts: true if returning as voltage
%   memm.type: returns the data type, as stored on disk
% By typing "disp(memm)" you can get a sense for other useful fields, some
% of which can be modified after creation.
%
% Various parameter/value pairs control the output:
%   tovolts (default value: true): if true, returns the data in units of
%     volts rather than A/D units;
%   blocksize (default value: 20Mb): allows you to specify the size of
%     the memory-mapped buffer, in bytes. Larger buffers consume more
%     memory but can yield faster performance.
%   contiguous (default false): if true, scan #s that are a 2-vector are
%     interpreted as a [min max] pair, and all values between them
%     (inclusive of the limits) them are returned.  This can yield higher
%     performance in cases where you want to read data in large blocks.
%     However, in applications like snippet cutting that may involve many
%     small blocks, you're better off with contiguous=false.
  
% Copyright 2006 by Timothy E. Holy

% Because need it earlier than other parameters, look for
% use_alt_header_file flag...
for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
        error(['After filename, input must consist of parameter name & ' ...
            'value pairs']);
    end
    switch lower(varargin{i})
        case 'use_alt_header_file'
            memm.use_alt_header_file = varargin{i+1};
        otherwise
            stoppoint = 1;
    end
end

  if isa(filename,'merecmm')
    % Copy constructor
    memm = filename;
  else
    memm.filename = filename;
    if isfield(memm,'use_alt_header_file')
        memm.header = readheader(memm.use_alt_header_file);
    else
        memm.header = readheader(filename);
    end
    % Set up lookup table for rapidly mapping channel #s to channel index
    maxchan = max(memm.header.channels);
    chan2ind = nan(1,maxchan+1);
    for i = 1:length(memm.header.channels)
      chan2ind(memm.header.channels(i)+1) = i; % +1 to account for chan0
    end
    memm.chan2ind = chan2ind;
    if (memm.header.minSample == 0)
      memm.type = 'uint16';
    else
      error('Fix me');
    end
    memm.offset = memm.header.headersize;
    memm.size = [memm.header.numch memm.header.nscans];
    [c,mxsz,endian] = computer;
    memm.swapbytes = ~(strcmp(lower(endian),lower(memm.header.endian)));
    memm.tovolts = true;
    memm.contiguous = false;
    blocksize = 20*1024^2;
    for i = 1:2:length(varargin)
      if ~ischar(varargin{i})
        error(['After filename, input must consist of parameter name & ' ...
               'value pairs']);
      end
      switch lower(varargin{i})
       case 'tovolts'
        memm.tovolts = varargin{i+1};
       case 'blocksize'
        blocksize = varargin{i+1};
       case 'contiguous'
        memm.contiguous = varargin{i+1};
       case 'use_alt_header_file'
        memm.use_alt_header_file = varargin{i+1};
       otherwise
        error(['Parameter ' varargin{i} ' not recognized.']);
      end
    end
    blocksz_in_scans = floor(blocksize/sizeof(memm.type)/ ...
                            memm.header.numch);
    if (blocksz_in_scans == 0)
      error('blocksize is too small for any values to be read');
    end
    memm.blocksz_in_scans = blocksz_in_scans;
    % Set up for conversion to volts; need this even if tovolts is false,
    % because the user might change tovolts later.
    memm.d2v_slope = (memm.header.voltageMax-memm.header.voltageMin)/ ...
        (memm.header.maxSample-memm.header.minSample);
    memm.d2v_offset = memm.header.voltageMax - ...
        memm.d2v_slope*memm.header.maxSample;
    % Finally, set up the memory mapping
    mm = memmapfile(filename,...
                    'format',{memm.type,...
                    [memm.header.numch min(memm.size(2),memm.blocksz_in_scans)],'x'},...
                    'repeat',1,'offset',memm.offset);
    memm.mm = mm;
    [memm.pathstr,memm.basestr,memm.extstr] = fileparts(mm.Filename);
    memm = class(memm,'merecmm');
  end