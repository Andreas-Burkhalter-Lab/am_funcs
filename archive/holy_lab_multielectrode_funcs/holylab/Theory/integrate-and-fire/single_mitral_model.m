function rate = single_mitral_model(I,t,mgp)
% Syntax:
%   rate = single_mitral_model(I,t,mgp)
% where
%   I is the injected current
%   t is the time window of the simultation
%   mgp is the mitral_granule parameters structure (see mitral_granule)

  p.n_cells = 1;
  p.dt = mgp.granule_timeconst/10; % take timesteps small compared to granule cell integration time
  p.conductance_fcn = @(v,t,tspike) mitral_granule(v,t,tspike,mgp);
  p.current_fcn = @(t) constant_current(t,I);
  [tspike,tsim,v] = integrate_fire(t,p);
  rate = length(tspike{1})/t;
  