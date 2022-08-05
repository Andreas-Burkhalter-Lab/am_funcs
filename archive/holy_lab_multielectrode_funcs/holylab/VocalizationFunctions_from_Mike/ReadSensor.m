function [transitions,initialVoltage,header] = ReadSensor(file,options)
% READSENSOR:  Read the proximity detection data from a file
%
% This function reads the header and proximity detection data from a file. 
% The read consists of the intial voltage followed by an index of when the 
% transitions from low to high or high to low occured. The read is not 
% processed in blocks.
%
% Syntax:
%   transitions = ReadSensor(file)
%   [transitions,initialVoltage,header] = ReadSensor(file)
%   snifintervals = ReadSensor(file,options)
% where
%    file is either a string and is treated as a filename or
%      is numeric and is treated as a file identifier
% and 
%    transitions is an index of when the transitions from low to high or 
%      high to low occured (by scan number)
%    initialVoltage is the voltage in the initial state (in volts)
%    header is the .det file header.
%
% If options are specified, you can instead get the output in the form of
% snifintervals, in which snifintervals is a 2-by-n matrix yielding the
% start and stop time of each period in which the mouse is poximal to the
% detector. options is a structure, with the following fields:
%   close: whether 'lo' or 'hi' indicates the nearness of the mouse,
%     depending on the detector. The capacitance detector should be 'lo',
%     the optical detector 'hi'. You must set this field.
%   intervals: set to true if you want these intervals (default true when 2
%     arguments are present)
%   tosecs: set to true if you want the times to be given in seconds rather
%     than scan intervals (default true with 2 arguments)
%   fixend: a parameter to handle a bug in certain .det files, in which the
%     number of scans was not properly adjusted for the infrequent sampling
%     of the proximity detection channel.  If present, nscans is divided by
%     fixend before any further processing. (default: 1)
%   
%
% See also: WRITESENSOR

if (nargin < 2)
  options.intervals = 0;
  options.tosecs = 0;
else
  if ~isstruct(options)
    error('options must be a structure');
  end
  if ~isfield(options,'intervals')
    options.intervals = 1;
  end
  if ~isfield(options,'tosecs')
    options.tosecs = 1;
  end
  if ~isfield(options,'close')
    error('You must supply a value for options.close');
  end
  if ~isfield(options,'fixend')
    options.fixend = 1;
  end
end

if (ischar(file))
  [fid,message] = fopen(file,'r');
  if (fid < 1)
    error(message);
  end
elseif (isnumeric(file))
  fid = file;
else
  error(['Do not recognize input ',file]);
end
% Read header
[header,headersize] = ReadAIHeaderWU1(fid);

initialVoltage = fread(fid,1,'float64');
transitions = fread(fid,header.numTransitions,'int32');


% Close file if "file" is a string
if (ischar(file))
  status = fclose(fid);
  if (status < 0)
    error('File did not close');
  end
end

% Parse into snifintervals matrix
if options.intervals
  t = [0;transitions;header.nscans/options.fixend];
  vthresh = 2.5;  % 2.5V is threshold for hi/lo
  if strncmp(options.close,'hi',2)
    v = (initialVoltage > vthresh);
  elseif strncmp(options.close,'lo',2)
    v = (initialVoltage <= vthresh);
  else
    error('options.close is neither lo nor hi');
  end
  v(2) = ~v;
  v = repmat(v,1,ceil(length(t)/2));
  %v(length(t)) = v(length(t)-1);  % The end point (last scan) is not a transition
  indx = find(v(1:length(t)-1)); % Since we need a stop time, we don't include the last scan
  transitions = [t(indx)';t(indx+1)'];
end

% Convert to secs, if needed
if options.tosecs
  transitions = transitions/header.scanrate;
end

  