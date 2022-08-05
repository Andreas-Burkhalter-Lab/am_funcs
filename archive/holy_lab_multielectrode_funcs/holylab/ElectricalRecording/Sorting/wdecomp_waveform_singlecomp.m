function w = wdecomp_waveform_singlecomp(v,t,trange,A)
% wdecomp_waveform_singlecomp: optimize a single filter given a signal, times, and amplitudes
% Syntax:
%   w = wdecomp_waveform_singlecomp(v,t,trange,A)
% where
%   v is a 1-by-n_samples signal
%   t is a 1-by-n_events set of event times (index to sample #)
%   trange is a 2-vector, [n_before n_after], giving the range of
%     waveform before and after each event that you seek to fit
%   A (optional, defaults to all 1) is a 1-by-n_events set of amplitudes
% and
%   w is a vector of length diff(trange)+1 that, when multiplied properly
%     by the amplitudes in A, best describes the signal v.
%
% See also: wdecomp_amplitude_rhs, wdecomp_amplitude_hessian.
 
% Copyright 2009 by Timothy E. Holy
  
  [n_channels,n_samples] = size(v);
  if (n_channels > 1)
    error('Multiple channels not yet implemented');
  end
  if (nargin < 4)
    A = ones(size(t));
  else
    if (size(A,1) > 1)
      error('Multiple waveform components not yet implemented');
    end
  end
  
  n_w = diff(trange)+1;
  
  % Trim t and A so that all spikes occur sufficiently away from the
  % edges
  keepFlag = t + trange(1) > 0 & t + trange(2) <= n_samples;
  t = t(keepFlag);
  A = A(keepFlag);

  % The equation for w will be Mw = b, where the matrix M has self-terms
  % and cross-terms (cross-terms comes when spikes overlap within trange)
  
  % Self-terms
  s = sum(A.^2);
  S = s*eye(n_w,n_w);
  
  % Cross-terms
  tmax = n_w;
  [tac,index] = autocorrspike(t,tmax);
  M = zeros(n_w,n_w);
  o = ones(1,n_w);
  for i = 1:length(tac)
    Atmp = A(index(1,i)) * A(index(2,i));
    Dtmp = diag(o(tac(i)+1:end),tac(i));
    M = M + Atmp * Dtmp;
  end
  M = M + M';
  M = M + S;
  
  % Right hand side
  b = zeros(1,n_w);
  rng = trange(1):trange(2);
  for i = 1:length(t)
    b = b + A(i)*v(t(i) + rng);
  end
  
  % Solve the linear equation
  w = (M\b')';
end
  
  
  