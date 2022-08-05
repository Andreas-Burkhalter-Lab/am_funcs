function [ephysout,utags] = ephystag(ephysin,tags,options)
% EPHYSTAG: assign tags to ephys structures
% Syntax:
%   ephysout = ephystag(ephysin,tags,options)
%   [ephysout,utags] = ephystag(....)
%
% where
%   ephysin is the input ephys structure array;
%   tags is a cell array of strings with the same number of elements as
%     ephysin, each string the name for the tag;
%   options is an optional structure array with the following fields:
%     addtime (default false): if true, will append the duration of valve
%       opening (in seconds) to the tag name, using a delimiter
%       timedelim;
%     timedelim (default $): the delimiter used before timing
%       information;
%     timeprec (default 0.025): if appending the time, first multiply the
%       duration by 1/timeprec, then round off, then multiply by
%       timeprec.  This insures that tiny differences (e.g., a few scans)
%       in open duration across trials do not end up receiving different
%       tags.
% and
%   ephysout is the output (tagged) ephys structure array (note it may
%     also have stimulus information loaded, if it wasn't already there);
%   utags (optional) is the list of unique tags, in alphabetical order,
%     useful for making sure that the tags assigned correspond with your
%     experimental design.
%
% See also: INTERVALSFROMSTIM, EPHYSSUBRANGE.
  
% Copyright 2004 Timothy E. Holy.
  
  if (nargin < 3)
    options = struct;
  end
  options = ephystagoptions(options);
  ephysout = ephysin;
  [ephysout.tag] = deal(tags{:});  % This is all that's needed for
                                   % many cases
  if options.addtime
    if ~isfield(ephysin,'stimulus')
      ephysout = ephysfetch(ephysout,'stimulus');
    end
    for i = 1:length(ephysout)
      nzindex = find(ephysout(i).stimulus(1,:)); % Find the first
                                                 % non-zero valve
      if (isempty(nzindex) | nzindex(1) == size(ephysout(i).stimulus,2))
        ontime = 0; % Default if you don't know what to do
      else
        ontime = diff(ephysout(i).stimulus(2,nzindex(1)+[0 1])) ...
                 / ephysout(i).scanrate;
        ontime = options.timeprec*round(ontime/options.timeprec); % Round
                                                                  % off times
      end
      ephysout(i).tag = [ephysout(i).tag,options.timedelim, ...
                         num2str(ontime)];
    end
  end
  if (nargout > 1)
    utags = unique({ephysout.tag});
  end
  
  
function options = ephystagoptions(options)
  if ~isfield(options,'addtime')
    options.addtime = 0;
  end
  if ~isfield(options,'timedelim')
    options.timedelim = '$';
  end
  if ~isfield(options,'timeprec')
    options.timeprec = 0.025;
  end
  