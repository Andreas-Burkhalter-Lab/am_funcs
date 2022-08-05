function y = sng2sound(sng,window)
% SNG2SOUND: re-create a waveform from the sonogram
% y = sng2sound(sng)
% sng may be computed by SPECGRAM or by SPARSESNG.
%
% See also: SPECGRAM, SPARSESNG, CLEANSOUND, FREQDIVIDE.

% See E. Moulines and J. Laroche, "Non-parametric techniques for
%   pitch-scale and time-scale modification of speech," Speech Communication
%   16: 175-205 (1995).
% Copyright 2001 by Timothy E. Holy <holy@pcg.wustl.edu>
  
  [nf,nt] = size(sng);
  nfft = 2*(nf-1);
  if (nargin < 2)
    window = hanning(nfft)';
  else 
    if (length(window) ~= nfft)
      error('window length and nfft do not match');
    end
  end
  npts = nfft/2*(nt+1);
  y = zeros(1,npts);
  tmp = zeros(1,nfft);
  for k = 1:nt
    cindx = 1+(k-1)*nfft/2:(k+1)*nfft/2;
    tmp(1:nf) = sng(:,k);
    tmp(nf+1:end) = conj(sng(end-1:-1:2,k));
    %if (k == 1)
    %        keyboard
    %end
    %mean(abs(imag(ifft(tmp))))
    if (k > 1 & k < nt)
      y(cindx) = y(cindx) + window.*real(ifft(tmp));
    else
      tmpwindow = window;
      if (k == 1)
        tmpwindow(1:nfft/2) = 1;
      end
      if (k == nt)
        tmpwindow(nfft/2:end) = 1;
      end
      y(cindx) = y(cindx) + tmpwindow.*real(ifft(tmp));
    end
  end
