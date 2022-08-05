function [r,tout] = rate1storder(t,stimulus,params)
% RATE1STORDER: compute firing rates from a 1st order binding model
% 
% Given receptor kinetics and stimulus concentrations, compute the
% occupancy fraction of the receptor.  The returned firing rate is
% proportional to the occupancy.
%
% More recently, this has been generalized to allow for spike rate
% adaptation.  The mean firing rate is computed over a window of past
% history, and the output firing rate is reduced multiplicatively, where
% the factor is linear in the previous firing rate history (after
% subtracting the spontaneous rate).
%
% Syntax:
%   r = rate1storder(t,stimulus,params)
% where
%   t is a vector of time points at which to evaluate the rate (in seconds;
%     must be sorted!)
%   stimulus is a 2-by-n matrix, stimulus(1,:) is the sequence of valve
%     identities (integers), and stimulus(2,:) is the sequence of
%     transition times (in seconds);
%   params is a structure with the following fields:
%     rmax: a vector (1-by-ncells) of peak firing rates (Hz)
%     rspont: a vector (1-by-ncells) of spontaneous firing rates (Hz)
%     cOverK: a matrix (nstim-by-ncells) giving c/K for each stimulus
%       (indexed by the valve number). Note that valve = 0 is set to give
%       cOverK = 0.
%     koff: a vector (1-by-ncells) of off rate constants (Hz)
%     inttime (optional): 1-by-ncells integration time constant for
%       measuring the spike rate over recent history, used only if you
%       want spike rate adaptation (note "integration" is done using
%       discrete time bins, so you'll get higher accuracy with more time
%       bins)
%     multfac (optional): 1-by-ncells fraction by which activity is
%       reduced if the time-averaged firing rate were rmax (multfac = 1
%       sets no adaptation)
% and
%   tout is a corrected times vector (to make sure all times fall within
%     the time window specified by the stimulus matrix);
%   r is a ncells-by-length(tout) matrix of firing rates.
%
% This function can be used to simulate spikes with EPHYSSIMULATEDSPIKES.
%
% See also: EPHYSSIMULATEDSPIKES, CONCCOMPARE.

ncells = length(params.rmax);
tindx = find(t >= stimulus(2,1) & t <= stimulus(2,end));
t = t(tindx);
tout = t;
for i = 1:ncells
  % First time period: assume stimulus has been in this state forever (at
  %   steady-state)
  tindx = find(t <= stimulus(2,2));
  if (stimulus(1,1) == 0)
    cK = 0;
  else
    cK = params.cOverK(stimulus(1,1),i);
  end
  pocc = cK/(1+cK);
  rlast = params.rmax(i)*pocc + params.rspont(i);
  r(i,tindx) = rlast;
  % Now handle transitions to new states
  for j = 2:(size(stimulus,2)-1)
    tstart = stimulus(2,j);
    tend = stimulus(2,j+1);
    tindx = find(t > tstart & t <= tend);
    if (stimulus(1,j) == 0)
      cK = 0;
    else
      cK = params.cOverK(stimulus(1,j),i);
    end
    kp = params.koff(i)*cK;
    ksum = kp+params.koff(i);
    expt = exp(-ksum*(t(tindx)-tstart));
    p = kp/ksum * (1 - expt) + pocc*expt;
    r(i,tindx) = params.rmax(i)*p + params.rspont(i);
    expt = exp(-ksum*(tend-tstart));
    pocc = kp/ksum*(1-expt) + pocc*expt;
  end
end
% If spike rate adaptation is desired, do it
%
% Note the adaptation factor is:
%   g(r) = (rmax - multfac*rspont - (1-multfac)*rbar)/(rmax-rspont)
% where rbar is the average rate over the integration time
%
if (isfield(params,'multfac') && isfield(params,'inttime'))
  for i = 1:ncells
    % As before, assume that the first period is at steady-state.
    tindx = find(t <= stimulus(2,2));
    if isempty(tindx)
      error(['Haven''t implemented spike rate adaptation if first epoch not ' ...
             'present.']);
    end
    r0 = r(i,1);
    rdiff = params.rmax(i)-params.rspont(i);
    rmdiff = params.rmax(i)-params.multfac(i)*params.rspont(i);
    rmean = r0*rmdiff/rdiff;
    rmean = rmean/(1 + (1-params.multfac(i))*r0/rdiff);
    r(i,tindx) = rmean;
    % Now handle transitions to new states
    updatefac = (tout(2)-tout(1))/params.inttime(i);
    for j = tindx(end)+1:size(r,2)
      rmean = (1-updatefac)*rmean + updatefac*r(i,j);
      r(i,j) = (rmdiff - (1-params.multfac(i))*rmean)/rdiff * r(i,j);
    end
  end
end