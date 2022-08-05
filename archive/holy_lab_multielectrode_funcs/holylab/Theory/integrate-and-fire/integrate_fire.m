function [tspike,t_out,v_out] = integrate_fire(trange,params)
% INTEGRATE_FIRE: compute spike times in arbitrary I-F networks
% An integrate-and-fire neuron can be described by the following
% differential equation:
%    I(t) = sum_j (V(t) - E_j)/R_j(t) + C dV/dt
% where I(t) is the time-varying injected current, V is the membrane
% voltage, R_j and E_j are the resistance and reversal potential of its
% conductances, and C is the capacitance. This equation is supplemented by
% the statement that the neuron fires when the membrane potential reaches
% Vth, at which point the membrane potential gets reset to 0.
%
% One can write this equation in a terms of units of time only by
% defining:
%     v = V/Vth (so the spiking threshold is 1)
%     e_j = E_j/Vth
%     i = I/qth, where qth = C Vth (has units of rate)
%     tau_j(t) = C R_j(t)
% yielding
%     dv/dt = i(t) - sum_j (v(t)-e_j)/tau_j(t).
% The neuronal model is specified by these dimensionally-reduced
% quantities.
%
% Implementation: note that we can write this equation as
%    dv/dt = iprime - v/tau,
% where
%   iprime = i + sum_j (e_j/tau_j)
%   1/tau = sum_j (1/tau_j)
% If tau and iprime are constants (not time-varying), then this equation
% can be solved exactly:
%   v = v0 + v1(1-exp(-t/tau)),
% where v1 = iprime*tau - v0. Note that the sign of v1 determines whether
% the membrane potential will rise or fall, and v0+v1 >= 1 implies that
% the neuron will eventually spike.
% Thus, the implementation allows the user to specify a parameter dt as a
% time interval over which one should assume tau and iprime are constant,
% and then solves the equations exactly over this time interval. The
% forward propagation halts, however, whenever any neuron spikes, since
% this could change the conductances.  This approach reduces numerical
% error compared with discretely propagating the equations forward in
% time, all the way to zero error in the (trivial) cases where the
% currents and conductances really are constant. Note, however, that the
% temporal step size will shrink as you have more an more neurons in your
% circuit, causing the simulation to take much longer.
%
% Syntax:
%   [tspike,t,v] = integrate_fire(trange,params)
% where
%   trange is either a scalar (the # of "seconds" to run the computation)
%     or a 2-vector [starttime stoptime], useful if you're wanting to
%     extend some previous results further in time;
%   params is a structure which may have the following fields:
%     n_cells (required): the # of cells to simulate
%     conductance_fcn (required): the function to call to determine the
%       tau_j(t) and e_j for each neuron. This function has syntax
%          conductance = conductance_fcn(v,t,tspike)
%       v is the current membrane voltage, t is the current time, and
%       tspike is a cell array containing the previous firing history of
%       the simulated neurons (so that the conductances may depend upon
%       firing). You have to write the actual conductance function
%       (examples are PASSIVE_NEURON_CONDUCTANCE and MITRAL_GRANULE). Any
%       parameters that need to be passed to the conductance function
%       should be supplied using an anonymous function.
%       The ouput, conductance, is a structure with two fields, tau (a
%       matrix, tau(i,j) is the jth time constant for the ith cell) and
%       reversal (reveral(i,j) is the jth reversal potential for the ith
%       cell).
%     current_fcn (required): the function to call to determine i(t) for
%       each neuron. This function has syntax
%         currents = current_fcn(t)
%       where t is the time. currents is a vector of size
%       n_cells-by-1. Again, any parameters needed by your current_fcn
%       should be passed by defining an anonymous function. An example is
%       CONSTANT_CURRENT.
%
%     v0 (default 0): the starting voltage for each of the cells (a
%       n_cells-by-1 vector)
%     tspike (default empty): any previous history of firing times (may be
%       needed for determining the conductances)
%     dt (default Inf): the time step to take forward. If this is set
%       to Inf, then the model will be integrated forward in time until
%       the next neuron spikes (or the end of the trange is reached),
%       assuming that all conductances and currents remain constant.
%       If you need either of these to change over time even in the
%       absence of spiking, then you should set dt to be something
%       relatively small on the scale of these changes.
%
% On output,
%   tspike is a cell array of spike times, one entry per neuron. If you are
%     temporally extending some previous results (by supplying a trange
%     and tspike option), those previous results will be concatenated with the
%     new results.
%   t is a vector of times at which the membrane voltages are evaluated
%     (this will _not_ be concatenated with previous history...);
%   v(i,j) is the membrane voltage of the ith neuron at time t(j).
%
% See also: PASSIVE_NEURON_CONDUCTANCE, MITRAL_GRANULE, CONSTANT_CURRENT.

% Copyright 2006 by Timothy E. Holy
  
  if (length(trange) < 2)
    trange(2) = trange(1);
    trange(1) = 0;
  end
  params = default(params,'v0',zeros(params.n_cells,1));
  params = default(params,'dt',Inf);
  params = default(params,'tspike',cell(params.n_cells,1));
  
  tspike = params.tspike;
  t = trange(1);
  v = params.v0;
  if (length(v) ~= params.n_cells)
    error('Mismatch between n_cells and v0');
  end
  t_out = zeros(1,0);
  v_out = zeros(params.n_cells,0);
  while (t < trange(2))
    % Calculate conductances & reversal potentials
    conductance = params.conductance_fcn(v,t,tspike);
    % Convert conductances & reversal potentials to effective currents
    % and time constants
    tau = 1./sum(1./conductance.tau,2);
    iprime = sum(conductance.reversal ./ conductance.tau, 2);
    % Add additional input currents
    iprime = iprime + params.current_fcn(t);
    % Determine the asymptotic voltage of each neuron, assuming tau and
    % iprime are constant
    v_asymp = iprime*tau; % this is v0+v1
    v1 = v_asymp - v;
    % Determine which (if any) neurons will eventually spike
    will_fire = (v_asymp > 1);
    % For those that will fire, determine the earliest spike time, and
    % step forward to that time (constrained by the maximum time step,
    % params.dt)
    t_fire = inf(params.n_cells,1);
    t_fire(will_fire) = tau(will_fire) .* log(v1(will_fire)./(v_asymp(will_fire)-1));
    dt = min(trange(2)-t,min(params.dt,t_fire));
    t = t+dt;
    if (nargout > 1)
      t_out(end+1) = t;
    end
    % Update the tspike register
    index_fired = find(t_fire == dt);  % handles "ties" (simultaneous firing)
    for i = 1:length(index_fired)
      tspike{index_fired(i)}(end+1) = t;
    end
    % Move the voltages forward to their values at the new time
    v = v - v1.*expm1(-dt./tau);
    v(index_fired) = 0;  % reset membrane voltage of firing cells
    if (nargout > 2)
      v_out(:,end+1) = v;
    end
  end
