function conductance = mitral_granule(v,t,tspike,params)
% MITRAL_GRANULE: inhibitory mitral conductances in a reciprocal network
% Ficticious granule cells temporally integrate the ouptut of mitral
% cells, then turn this activity back onto the mitral cells as
% (instantaneous) shunting conductances.
%
% Syntax:
%   conductance = mitral_granule(v,t,tspike,params)
% where
%   v,t,tspike are all described in INTEGRATE_FIRE;
%   params has the following fields:
%     granule_timeconst: the time constant over which spiking is integrated
%       in "granule cells"
%     inhib_matrix: the auto-inhibitory matrix, i.e., inhib_matrix(i,j)
%       describes the effective strength of the jth mitral cell in
%       driving inhibition for the ith mitral cell (in units of
%       1/tau/firing rate, which means it's dimensionless)
%     external_tau (optional): a constant (or column vector) of
%       additional inhibitory taus; leak conductances can be included here
%     min_tau (optional): a minimum tau, allowing for granule cell output
%       to saturate.
%
% See also: INTEGRATE_FIRE.

% Copyright Timothy E. Holy
  
  % Integrate the history of mitral cell firing; use cellfun to try to
  % speed this up
  if (length(v) > 1)
    filter_func = @(x) mg_filtspike(x,t,params.granule_timeconst);
    filtspike = cellfun(filter_func,tspike);
  else
    filtspike = mg_filtspike(tspike{1},t,params.granule_timeconst);
  end
  % Apply this integrated activity as inhibitory conductances
  inv_tau = params.inhib_matrix * filtspike;
  % Add any extra shunting conductances
  if isfield(params,'external_tau')
    inv_tau = inv_tau + 1./params.external_tau;
  end
  conductance.tau = 1./inv_tau;
  if isfield(params,'min_tau')
    % Saturation of the inhibitory input
%     min_tau = params.min_tau;
%     if isscalar(min_tau)
%       min_tau = min_tau * ones(size(conductance.tau));
%     end
    %conductance.tau = max(conductance.tau,params.min_tau);  % hard limit
    %conductance.tau = min_tau .* sqrt(1+(conductance.tau ./ min_tau).^2); % old soft limit
    conductance.tau = conductance.tau + params.min_tau;
  end
  conductance.reversal = 0;
  
function filtspike = mg_filtspike(tspike,t,timeconst)
% Filter the spike history with an exponential function
  t_earliest = t-10*timeconst;
  index = find(tspike < t_earliest,1,'last');
  if isempty(index)
    index = 0;
  end
  dt = t-tspike(index+1:end);
  filtspike = sum(exp(-dt/timeconst));
  