function wave = template_copies(template,t,a,trange)
% TEMPLATE_COPIES: place variable-amplitude copies of a template at defined times
%
% Syntax:
%   wave = template_copies(template,t,a,trange)
% where
%   template is a vector of values of a "base waveform"
%   t is a 1-by-n vector of "times" at which to place the base waveform
%     (in index units, although fractional times are supported)
%   a is a 1-by-n vector containing the amplitude to apply for each event
%   trange is a 2-vector, [tmin tmax], over which to generate the output
% and
%   wave is a vector containing the output waveform.

% Copyright 2009 by Timothy E. Holy

  n_events = length(t);
  wave = zeros(1,diff(trange)+1);
  rng0 = 0:length(template)-1;
  if (size(template,1) > 1)
    wave = wave';
    rng0 = rng0';
  end
  for i = 1:n_events
    % Calculate the integer and fractional part of the offset
    ti = round(t(i));
    tf = t(i)-ti;
    % Determine which indices lie within the range of I
    x = ti+rng0;
    keepFlag = x >= trange(1) & x <= trange(2);
    template_shift = interp1(rng0,template,rng0(keepFlag)-tf);
    template_shift(isnan(template_shift)) = 0;
    xo = x(keepFlag)-trange(1)+1;
    wave(xo) = wave(xo) + a(i)*template_shift;
  end
