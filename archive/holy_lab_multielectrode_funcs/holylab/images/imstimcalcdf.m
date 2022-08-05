function ipo = imstimcalcdf(ip,options)
% IMSTIMCALCDF: image fluorescence changes upon stimulation
%
% This function looks for valve transitions, and then selects a few frame
% before the transition to be the background, a few after the transition
% to be the foreground.  By default, it calculates 
%           deltaf/f  =  (foreground-background)/background
% and returns the results.
%
% Syntax:
%   ipout = imstimcalcdf(ip,options)
% where
%   ip is an imphys structure;
%   options is a structure with the following fields:
%     flush (default 0): the number of the flush valve;
%     numbg (default 4): the number of frames used in computing the
%       background (immediately before the stimulus turns on);
%     numfg (default 4): the number of frames used in computing the
%       foreground (backward from the stimulus offset);
%     valvenum (default all): valves to analyze;
%     nonorm: if true, doesn't normalize the images (calculates deltaf
%       rather than deltaf/f);
%     filter: if present, applied on the average of foreground, and the
%             average of background
%     shiftfg, shiftbg: the #frames to shift when count fg and bg
% and
%   ipout is an imphys structure containing the results.
%
% See also: IMSTIMVIEWDF.

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
    
  % Find the stimulus onsets & offsets
  stim = [ip.stimulus];
  indx_on = find(stim(1:end-1) == options.flush & ...
                 stim(2:end) ~= options.flush) + 1;
  indx_off = find(stim(1:end-1) ~= options.flush & ...
                  stim(2:end) == options.flush) + 1;
  if (stim(end) ~= options.flush)
    indx_off = [indx_off length(stim)+1];
  end
  vlv_on = stim(indx_on);
  indxkeep = indexainb(vlv_on,valvenum);
  indx_on = indx_on(indxkeep);
  indx_off = indx_off(indxkeep);
  if (length(indx_on) ~= length(indx_off))
    'Valvenumbers:'
    stim
    error('Number of onsets is not equal to the number of offsets');
  end
  
  % Loop over transitions
  ntrans = length(indx_on);
  for i = 1:ntrans
    % Fill in stimulus & timing info
    tmp = ip(indx_on(i));
    % Compute averaged & background-subtracted image
    rng_bg = indx_on(i)-options.numbg:indx_on(i)-1;
    rng_bg=rng_bg+options.shiftbg;
    rng_fg = indx_off(i)-options.numfg:indx_off(i)-1;
    rng_fg=rng_fg+options.shiftfg;
    imbg = {ip(rng_bg).image};
    imfg = {ip(rng_fg).image};
    bg = immean(imbg{:});
    fg = immean(imfg{:});
    
    % filter both foreground and background before do division:
    if(isfield(options,'filter'))
       bg=imfilter(bg, options.filter);
       fg=imfilter(fg, options.filter);
    end
    
    if options.nonorm
      tmp.image = fg - bg;
    else
      tmp.image = (fg-bg)./bg;
    end
    % Record the stack numbers
    tmp.stacknumbg = [ip(rng_bg).stacknum];
    tmp.stacknumfg = [ip(rng_fg).stacknum];
    
    % the name of computation:
    if(options.nonorm)
       tmp.computation='\DeltaF';
    else
       tmp.computation='\DeltaF/F';
    end
    ipo(i) = tmp;
  end
  
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
  if ~isfield(options,'nonorm')
    options.nonorm = 0;
  end
  if(~isfield(options, 'shiftfg'))
     options.shiftfg=0;
  end
  if(~isfield(options, 'shiftbg'))
     options.shiftbg=0;
  end
  