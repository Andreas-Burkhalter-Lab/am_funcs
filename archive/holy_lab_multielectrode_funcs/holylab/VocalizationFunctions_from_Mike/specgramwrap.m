function [B,xoldo] = specgramwrap(x,xold,nfft,window)
% SPECGRAMWRAP: specgram for long signals
%
% This function is useful for computing the spectrogram of signals that
% may not be held in memory all at once.  The signal may be processed in
% blocks; the boundaries between blocks are handled properly by the xold
% arguments to this function.
%
% Syntax:
%   [B,xoldo] = specgramwrap(x,xold,nfft,window)
% where
%   x is the current signal, with a number of points equal to an integer
%     multiple of nfft/2;
%   xold contains past data, with length nfft/2.  On the first call, the
%     user should set this to zeros, for zero-padding.  On all other
%     calls, the user should pass the output xoldo from the previous
%     call to this function.
%   nfft (optional) is the number of points in the fourier transform;
%     this must be a power of 2
%     (default value: 256)
%   window (optional) is the windowing function (default: hanning(nfft))
% and
%   B is the specgram matrix, B(freq,time).  It has size
%    (nfft/2+1)-by-(2*(nsamp-nfft/2)/nfft).
%   xoldo contains the last nfft/2 points of x
%
% The overlap between successive columns of B is nfft/2.  Real signals
%   are assumed.
% This function gives the same output as
%   specgram([xold x], nfft).
% However, this version is somewhat faster than specgram.
%
% The only subtlety concerns the proper handling of the tails of the
% signal; see INVSPECGRAMWRAP for information on the proper way to
% analyze the data to insure that all data can be reconstructed.
%
% See also: SPECGRAM, INVSPECGRAMWRAP, SPARSESNG.

% Copyright 2002 by Timothy E. Holy <holy@pcg.wustl.edu>

% Fill in default arguments and do error checking on inputs
  if (nargin < 3)
    nfft = 256;
  end
  if (nargin < 4)
    window = hanning(nfft);
  end
  l2 = log2(nfft);
  if (floor(l2) ~= l2)
    error('nfft must be a power of 2');
  end
  ncol = length(x)/nfft * 2;
  if (floor(ncol) ~= ncol)
    error('x must have a length which is an integer multiple of nfft/2');
  end
  if (length(xold) ~= nfft/2)
    error('xold must have length nfft/2');
  end
  
  % Reshape the data to make columns of size nfft,
  % with proper overlap and old data at the beginning
  % This can be accomplished simply by shifting and stacking copies of
  % the reshaped waveform
  B = reshape(x,nfft/2,ncol);
  B = [xold(:) B(:,1:end-1); B];
  
  % Window the data
  B = B .* repmat(window,1,ncol);

  % Fourier transform
  B = fft(B,nfft);
  
  % Assume input signal is real, and take only the necessary part
  B = B(1:nfft/2+1,:);

  % Prepare for the next call
  xoldo = x(end-nfft/2+1:end);
