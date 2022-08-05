function w = wdecomp_waveform(v,t,trange,A)
% wdecomp_waveform: optimize kernels given a signal, times, and amplitudes
% Syntax:
%   w = wdecomp_waveform(v,t,trange,A)
% where
%   v is an n_channels-by-n_samples signal
%   t is a 1-by-n_events set of event times (an index of sample #s)
%   trange is a 2-vector, [n_before n_after], giving the range of
%     waveform before and after each event that you seek to fit
%   A is a n_components-by-n_events set of amplitudes. If you omit this
%     it will default to a single component with the amplitude of each
%     event set to 1.
% Then,
%   w is a n_channels-by-len-by-n_components array (len =
%     diff(trange)+1), for which w(:,:,i) is the ith component.  When
%     multiplied by the amplitudes in A, and placed at the times
%     indicated by t, this best describes the signal v.
%
% See also: wdecomp_amplitude_rhs, wdecomp_amplitude_hessian.
 
% Copyright 2009 by Timothy E. Holy
  
  [n_channels,n_samples] = size(v);
  if (nargin < 4)
    A = ones(size(t));
  end
  [n_components,n_events] = size(A);
  if (length(t) ~= n_events)
    error('Mismatch between t and A');
  end
  
  n_w = diff(trange)+1;
  
  % Trim t and A so that all spikes occur sufficiently away from the
  % edges
  keepFlag = t + trange(1) > 0 & t + trange(2) <= n_samples;
  t = t(keepFlag);
  A = A(:,keepFlag);

  % The equation for w will be Mw = b
  
  % Right hand side

  b = zeros(n_channels,n_w,n_components);
  rng = trange(1):trange(2);
  for tIndex = 1:length(t)
    vsnip = v(:,t(tIndex) + rng);
    for cIndex = 1:n_components
      b(:,:,cIndex) = b(:,:,cIndex) + A(cIndex,tIndex)*vsnip;
    end
  end
  
  % There are no cross-terms on the channels, so permute b so that the
  % matrix will be block-diagonal
  bperm = permute(b,[3 2 1]);  % now in comp-time-channel order
  
  % The left hand side
  % M has "self terms" (spikes interacting with self) and "cross terms"
  % (spikes interacting with other spikes)

  % Self-term coefficients
  S_coefs = A*A';
  
  % Cross-term coefficients
  tmax = n_w;
  [tac,index] = autocorrspike(t,tmax);
  cl_tac = agglabel(tac);  % collect all event-pairs that have a given time difference
  n_cl = length(cl_tac);
  C_coefs = zeros(n_components,n_components,n_cl);
  for i = 1:n_cl
    this_timegap = cl_tac{i};
    C_coefs(:,:,i) = A(:,index(1,this_timegap)) * ...
	A(:,index(2,this_timegap))';
  end
  
  % Form the cross-term matrix
  n_tot = n_w*n_components;
  C = zeros(n_tot,n_tot);
  for i = 1:n_cl
    C = C + wdcw_place_on_diag(C_coefs(:,:,i),i,n_w);
  end
  C = C + C';
  
  % Form the self-term matrix
  S = wdcw_place_on_diag(S_coefs,0,n_w);
  
  M = S + C;
  
  % Solve the linear equation & put in expected form
  w = M\reshape(bperm,[n_components*n_w n_channels]);
  w = reshape(w,[n_components n_w n_channels]);
  w = ipermute(w,[3 2 1]);
end
  
  
function D = wdcw_place_on_diag(C,i,nt)
  Csz = size(C);
  Crep = repmat({C},1,nt-i);
  D = blkdiag(Crep{:});
  Dz = zeros(Csz(1)*(nt-i),Csz(2)*i);
  D = [Dz D; zeros(Csz*i) Dz'];
end
