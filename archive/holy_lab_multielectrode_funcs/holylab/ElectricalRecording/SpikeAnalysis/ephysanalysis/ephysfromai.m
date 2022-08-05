function ephysout = ephysfromai(filenames,options)
% EPHYSFROMAI: start an ephys structure from raw data files
%
% The header of the files is parsed to set the initial data
% in the ephys structure.  No "real" data is loaded; the structure is
% simply prepared for future manipulation.
%
% Syntax:
% ephysout = ephysfromai(filenames)
% ephysout = ephysfromai(headers)
%
% ephysout = ephysfromai(...,options)
%
% where
%   filenames is a string (or cell array of strings) with the filename
%     of the waveform, envelope, or snippet file;
%   headers is the header (or structure array of headers) read by
%     READHEADER.
%   options is an optional structure argument:
%     options.machfmt controls the opening format of the file
%     options.scanvalves, when true, causes the userheader to be parsed to
%       try to assign a label to each valve (default true)
%     options.autofn, when true, searches the current directory for files
%       with the same basename and conventional extensions (.merec, .env,
%       .vlv, and .ssnp) to fill in the filenames (default true)
%     options.usefullpath, when true, switches from searching the current
%       directory for filenames to fill in to the original path
%
% See also: READHEADER.
  
% Changes: 
%   TEH, 05-03-2004
%     Eliminate the "wholeheader" field from the header (to save memory)
%     Add the auto filename option
%
%   RCH, 08-15-2005
%     Modified auto filename system so that if filenames are of the 
%     form "*_clean.ssnp", "_clean" will be dropped from the basename
%   RCH, 10-25-2005
%     Strengthened the dropping of _clean from everything but .ssnp
%     filename
%   RCH, 10-20-2006
%     Added ability to recognize when it's being asked to process a
%     .ECG_ssnp 'snippet' file and enter the snipfile name appropriately
%   RCH, 02-11-2007
%     Added usefullpath option

  if (nargin < 2 | ~isfield(options,'scanvalves'))
    options.scanvalves = 1;
  end
  if ~isfield(options,'autofn')
    options.autofn = 1;
  end
  options = default(options,'usefullpath',true);
  options = default(options,'lfp_flag',0);

  % Handle the two modes of calling: headers or filenames
  havenames = 0;
  if (ischar(filenames) | (iscell(filenames) & ischar(filenames{1})))
    if ischar(filenames)
      filenames = {filenames};
    end
    nfiles = length(filenames);
    for i = 1:nfiles
      h(i) = readheader(filenames{i},options);  % Read in the headers
    end
    havenames = 1;
  elseif isstruct(filename)
    h = filename;
    nfiles = length(h);
  else
    error('Illegal argument');
  end

  % Fill in the basic elements
  [ephysout(1:nfiles).header] = deal([]);   % Creates the structure with blank field
  [ephysout.cellnums] = deal('allonchan');
  [ephysout.channels] = deal(h.channels);
  [ephysout.scanrate] = deal(h.scanrate);
  [ephysout.tovolts] = deal(h.scalemult);
  
  % Handle endian issues (default copies to all)
  if isfield(options,'machfmt')
    [ephysout.wavefilemachfmt] = deal(options.machfmt);
    [ephysout.envelopefilemachfmt] = deal(options.machfmt);
    [ephysout.snipfilemachfmt] = deal(options.machfmt);
  end
  
  % Do any necessary parsing
  for i = 1:nfiles
    if isfield(h(i),'wholeheader')
      htemp = rmfield(h(i),'wholeheader');
      ephysout(i).header = htemp;
    else
      ephysout(i).header = h(i);
    end
    if (isfield(options,'scanvalves') & options.scanvalves)
      [txt,vlabel] = ParseUsrHdr(h(i).usrhdr);
      ephysout(i).valvelabels = vlabel;
    end
    ephysout(i).scanrange = [1 h(i).nscans];
    if havenames
      [pth,basename,ext] = fileparts(filenames{i});
%       if isempty(pth)
%           pth = pwd;
%       end
      ephysout(i).basefilename = basename;
      if options.autofn
        if length(strfind(basename,'_clean')) > 0
          lengthBasenameOriginal = strfind(basename,'_clean')-1;
          basenameSnips = basename;
          basenameOriginal = basename(1:lengthBasenameOriginal);
        else
          basenameSnips = basename;
          basenameOriginal = basename;
        end
        if (options.usefullpath && ~isempty(pth))
            base = [pth filesep basenameOriginal];
            baseSnips = [pth filesep basenameSnips];
        else
            base = basenameOriginal;
            baseSnips = basenameSnips;
        end
        thefileinfo = dir([base '.vlv']);
        if ~isempty(thefileinfo)
          ephysout(i).stimulusfile = [base '.vlv'];
        end
        thefileinfo = dir([base '.merec']);
        if ~isempty(thefileinfo)
          ephysout(i).wavefile = [base '.merec'];
        end
        thefileinfo = dir([base '.env']);
        if ~isempty(thefileinfo)
          ephysout(i).envelopefile = [base '.env'];
        end
        if strcmp(ext,'.ssnp')
          thefileinfo = dir([baseSnips '.ssnp']);
          if ~isempty(thefileinfo)
            ephysout(i).snipfile = [baseSnips '.ssnp'];
          end
        elseif options.lfp_flag
            ephysout(i).snipfile = [baseSnips '.ssnp'];
                % because will actually want to get real snipinfo
        else
            ephysout(i).snipfile = [baseSnips '.ECG_ssnp'];
        end
        ephysout(i).basefilename = basenameOriginal;
      end
    end
  end
    