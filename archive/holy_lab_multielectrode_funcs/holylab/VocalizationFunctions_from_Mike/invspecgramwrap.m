function [x,xhalfo] = invspecgramwrap(sng,xhalf,window)
% INVSPECGRAMWRAP: re-create a waveform from the sonogram
% This function permits reconstruction even when the signal/sonogram is
% too large to fit in memory simultaneously, requiring processing in
% chunks.
%
% Syntax:
%   [x,xhalfo] = invspecgramwrap(sng,xhalf,window)
% where
%   sng is computed by SPECGRAMWRAP, and may be sparse.
%   xhalf is the previous tail of the x output, with length nfft/2
%     (initial to zeros, then on successive calls pass xhalfo as the input)
%   window (optional) is the window function (default: hanning(nfft))
% and
%   x is the reconstructed signal for this block
%   xhalfo is "half" of the signal for the first nfft/2 points of the
%     next block
% (This assumes half-overlap.)
%  
% The tails of the signal can be truncated or corrupted unless you treat
% these specially.  The following example shows how you can be sure to
% capture these tails properly:
%   Computing the spectrogram: suppose blocks are of the size n*nfft/2,
%   and define zer = zeros(1,nfft/2).  Then the following code will save
%   all the relevant data:
%       xold = zer;
%       for i = 1:nblocks
%         [B,xold] = specgramwrap(x,xold,nfft);
%         (Now write B to file or concatenate in memory)
%       end
%       % Do one extra column
%       [B1,xold] = specgramwrap(zer,xold,nfft);
%       (Append this B1 to the file/memory)
%
%   Inverting the spectrogram:
%       (B1 = first _column_ of the specgram file)
%       [x,xhalf] = invspecgram(B1,zer);
%       (You can throw away x, our purpose was to prepare xhalf)
%       for i = 1:nblocks
%         (B = next block from file---don't re-read the first column!)
%         [x,xhalf] = invspecgram(B,xhalf);
%         (x is now the signal for this block; you could save to disk or
%            concatenate in memory) 
%       end
%
% See also: SPECGRAMWRAP, SPECGRAM, SPARSESNG, CLEANSOUND, FREQDIVIDE.

% See E. Moulines and J. Laroche, "Non-parametric techniques for
%   pitch-scale and time-scale modification of speech," Speech Communication
%   16: 175-205 (1995).
% Copyright 5-03-01 by Tim Holy
% Changelog 02-21-02 (Tim Holy)
%   Made it work with processing in chunks.
  [nf,nt] = size(sng);
  nfft = 2*(nf-1);
  if (nargin < 3)
    window = hanning(nfft);
  else 
    if (length(window) ~= nfft)
      error('window length and nfft do not match');
    end
  end
  % Compute the window normalization factor
  wnorm = window(1:nfft/2) + window(nfft/2+1:nfft);
  % Re-form the complete spectrogram
  % (we assumed real signals, so we could throw away
  %  half of the spectrogram)
  x = [sng; conj(sng(end-1:-1:2,:))];
  % Inverse fourier transform
  x = ifft(x);
  % Preserve the tail for the next iteration
  xhalfo = x(nfft/2+1:nfft,end);
  % Now re-synthesize by overlap-add
  x = x(1:nfft/2,:) + [xhalf, x(nfft/2+1:nfft,1:end-1)];
  % Fix the normalization
  x = x ./ repmat(wnorm,1,nt);
  % Now turn back into a vector
  x = real(reshape(x,1,prod(size(x))));
