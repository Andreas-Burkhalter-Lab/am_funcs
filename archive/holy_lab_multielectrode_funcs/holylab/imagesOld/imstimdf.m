function ipp = imstimdf(ip,options)
% IMSTIMDF: image fluorescence changes upon stimulation
% Syntax:
%   ipp = imstimdf(ip,options)
% where
%   ip is an imphys structure;
%   options is a structure with the following fields:
%     flush (default 0): the number of the flush valve;
%     numbg (default ?): the number of frames used in computing the
%       background (immediately before the stimulus turns on);
%     numfg (default 4): the number of frames used in computing the
%       foreground (backward from the stimulus offset);
%     valvenum (default all): valves to analyze;
%     spatialwidth (default 8): the width of the averaging filter in
%       pixels; if this is empty, no spatial filtering is performed;
%     nogui: if set, operates in "batch" mode;
% and
%   ipp is an imphysp structure (see IMPHYSP) containing the results.
%
% See also: IMPHYSP, IMPHYSPROCESS, IMSHOWDF.

% Copyright 2005 by Timothy E. Holy & Jason Guo

  % Parse options
  if (nargin < 2)
    options = struct;
  end
  options = ispoptions(ip,options);
  if isfield(options,'valvenum')
    valvenum = options.valvenum;
  else
    valvenum = setdiff([ip.stimulus],options.flush);
  end
    
  % Show an image to get crop region
  if ~isfield(options,'nogui')
    hCropRect = imselrect(imphysfetch(ip(1)));
    ip = imcroprect(ip,hCropRect);
  end
  
  % Find the stimulus onsets & offsets
  stim = [ip.stimulus];
  indx_on = find(stim(1:end-1) == options.flush & ...
                 stim(2:end) ~= options.flush);
  indx_off = find(stim(1:end-1) ~= options.flush & ...
                  stim(2:end) == options.flush);
  if (length(indx_on) ~= length(indx_off))
    'Valvenumbers:'
    stim
    error('Number of onsets is not equal to the number of offsets');
  end
  
  % Loop over transitions
  ntrans = length(indx_on);
  for i = 1:ntrans
    rng_bg = indx_on(i)-options.numbg+1:indx_on(i);
    rng_fg = indx_off(i)-options.numfg+1:indx_off(i);
    ipp(i).oimfile = imphyscopy(ip([rng_bg rng_fg]),'imfile','check');
    ipp(i).ostacknum = [ip([rng_bg rng_fg]).stacknum];
    ipp(i).ostackweight = [repmat(-1/options.numbg,1,options.numbg) ...
                        repmat(1/options.numfg,1,options.numfg)];
    %ipp(i).filterwidth = options.spatialwidth;
    %ipp(i).resampmag = 0.5;
  end
  
  ipp = imphysprocess(ip,ipp);
  
  
function options = ispoptions(ip,options)
  if ~isfield(options,'flush')
    options.flush = 0;
  end
  if ~isfield(options,'numbg')
    options.numbg = 4;
  end
  if ~isfield(options,'numfg')
    options.numfg = 4;
  end
  if ~isfield(options,'spatialwidth')
    options.spatialwidth = 8;
  end