function t = spikesim(params,N)
% SPIKESIM: create simulated spikes for a refractory bursting neuron
%
% The neuron's spike-triggered firing rate is of the following form:
%   r = 0    if    0 <= t <= t1;
%   r = r1   if    t1 < t <= t2;
%   r = r2   if    t2 < t.
%
% Syntax:
%   t = spikesim(params,N)
% where
%   params is a structure with the fields r1,r2,t1,t2;
%   N is the desired number of simulated spikes
% and
%   t is a vector of interspike intervals.
%
% An alternative syntax is
%   meantime = spikesim(params),
% which returns the (scalar) mean interspike interval for the firing rate
% model specified by params.
  
  % See if we only want to return the mean time
  if (nargin < 2)
    pnoburst = exp(-params.r1*(params.t2-params.t1));
    t = params.t1 + (1-pnoburst)/params.r1 + pnoburst/params.r2;
    return
  end
  % OK, generate a number of random interspike intervals
  s = rand(1,N);
  lns = -log(s);
  Rcut = params.r1*(params.t2-params.t1);
  indx = find(lns < Rcut);
  t = zeros(1,N);
  t(indx) = lns(indx)/params.r1 + params.t1;
  indx = find(lns >= Rcut);
  t(indx) = (lns(indx) - Rcut)/params.r2 + params.t2;

