function [rmax,tmax,rmax_pre,tmax_pre] = deltarmax_calculation(t,options)
% DELTARMAX_CALCULATION: computing maximum firing rate in a response
% Syntax:
%   [rmax,tmax] = deltarmax_calculation(t,options)
%   [rmax,tmax,rmax_pre,tmax_pre] = deltarmax_calculation(t,options)
% where
%   t is a vector of spike times.  The units of t specify the units of
%     all other temporal quantities (i.e., if t is measured in seconds,
%     then all times are in seconds and all rates are in Hz)
%   options contains the fields
%     tmin in seconds (default 0). is the minimum amount of time over which a firing
%       rate must be measured
%     pre in seconds. This must be specified if you want to calculate the
%       rmax during the pre-stimulus interval (rmax_pre). For example, if
%       pre = -10, will look over the 10 sec prior to stim delvery to
%       calculate rmax_pre. 
% and
%   rmax is the maximum average firing rate measured from 0.  That is, it
%     is the number of spikes occuring before time tmax, divided by tmax.
%   tmax is the "cutoff" time for which the maximum average firing rate
%     was calculated.
%  
%11/2008 HAA added options field to include tmin and added pre.  


% See also: DELTARMAX.
  
% Copyright 2008 by Timothy E. Holy
  
  if (nargin < 2)
  options = struct;
  end
  
  options = default(options,'tmin',1);
  tmin = options.tmin;
   
 
  t_pre = t(t<0); % seperate out spikes that arrive before stimulus delivery
  t = t(t > 0); % discard spikes that arrived before time 0
  
  [rmax,tmax] = drmc_rmax(t,tmin);
  
  if (nargout > 2)
    % Check to see that the user specified "pre")
    if ~isfield(options,'pre')
      error('When calculating pre-stimulus firing rate, need to specify the starting time');
    end
    [rmax_pre,tmax_pre] = drmc_rmax(t_pre-options.pre,tmin);
  end

  
function [rmax,tmax] = drmc_rmax(t,tmin)  
  if isempty(t)
    % There are no spikes!
    rmax = 0;
    tmax = tmin;
    return
  end
  
  n = 1:length(t);  % calculate the cumulative # of spikes
  r = n(:)./t(:);         % calculate the average rate at the time of each spike
  % Choose the spiketime at a time bigger than tmin that yields the
  % largest firing rate
  [rmax,maxIndex] = max(r(t > tmin));
  tmax = t(maxIndex + sum(t <= tmin));
  if (tmin > 0)
    % It's possible that the largest average firing rate would be
    % measured at the tmin cutoff; compare what we have to what we'd get
    % there
    rmax_at_tmin = sum(t <= tmin)/tmin;
    if (isempty(rmax) || rmax < rmax_at_tmin)
      rmax = rmax_at_tmin;
      tmax = tmin;
    end
  end
    
  
