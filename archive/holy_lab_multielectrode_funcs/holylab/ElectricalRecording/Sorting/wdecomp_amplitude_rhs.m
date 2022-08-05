function b = wdecomp_amplitude_rhs(v,w,t)
% WDECOMP_AMPLITUDE_RHS: waveform decomposition for the amplitudes
% Computes the projections of the voltage waveform onto the components.
% Syntax:
%   b = wdecomp_amplitude_rhs(v,w,t)
% where
%   v is an n_channels-by-n_samples matrix containing the voltage
%     waveform;
%   w is an n_channles-by-l-by-n_components array, where w(i,j,k) is the
%     value of the kth component at time offet j on channel
%     i. Alternatively, if you're processing data in blocks you can
%     achieve faster performance (at the cost of increased memory usage)
%     by instead supplying conj(wfft), the conjugate of the fourier
%     transform of a suitably-padded w (must have l = n_samples in that
%     case).  Because each component will therefore be of the size of the
%     data v, you may need to process data in relatively small blocks.
%   t is a vector of spike times (as an index, i.e., scan number),
%     indicating the time of the beginning (_not_ peak time) of each
%     event. In terms of the peak time pt, t = pt - lleft, where lleft is
%     the number of samples in the components to the left of the peaks.
% and
%   b is the dot product of the voltage waveform with the components, an
%     n_components-by-n_times matrix.
%
% See also: WDECOMP_AMPLITUDE_HESSIAN.
  
% Copyright 2007 by Timothy E. Holy
  
  [n_channels,l,n_components] = size(w);
  if (size(v,1) ~= n_channels)
    error(['The voltage and the components do not have the same number of' ...
	   ' channels']);
  end
  n_scans = size(v,2);
  n_times = length(t);
  % Compute the correlations using the Fourier transform
  % That way, we have guaranteed performance no matter how many times
  % in t there are.
  vfft = fft(v,[],2);
  b = zeros(n_components,n_times);
  % Handle the components one at a time, so that we don't overflow memory
  for compIndex = 1:n_components
    wtmp = w(:,:,compIndex);
    if isreal(wtmp)
      % Components were supplied in temporal form
      if (l < n_scans)
        % Pad component with zeros
        wtmp(:,n_scans) = 0;
      end
      wfft = conj(fft(wtmp,[],2));  % Take the conjugate for correlation
    else
      % Components were already supplied as conj(fft(w))
      wfft = wtmp;
    end
    pfft = sum(wfft.*vfft,1);
    p = ifft(pfft);
    b(compIndex,:) = p(t);
  end
  