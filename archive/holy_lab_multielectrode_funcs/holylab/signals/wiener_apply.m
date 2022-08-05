function u = wiener_apply(v,w)
% WIENER_APPLY: apply a Wiener filter to a signal
% Syntax:
%   u = wiener_apply(v,w)
% where
%   v is the signal (of size n_samples-by-n_channels)
%   w is the Wiener filter(s) (of size
%          n_samples-by-n_channels-by-n_templates)
% and
%   u is the output (of size n_samples-by-1-by-n_templates)
%
% See also: WIENER.

% Copyright 2007 by Timothy E. Holy

  if (size(v,1) < size(w,1))
    v(size(w,1),:) = 0;
  else
    w = wiener_resize(w,size(v,1));
  end
  wfft = fft(w);
  vfft = fft(v);
  u = zeros([size(v,1) 1 size(w,3)],class(v));
  for templateIndex = 1:size(w,3)
    ufft = sum(wfft(:,:,templateIndex) .* vfft,2);
    u(:,templateIndex) = ifft(ufft);
  end
  