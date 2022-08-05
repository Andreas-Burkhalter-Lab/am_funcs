function fftinfo = register_fftprep(im_sz,sigma)
% REGISTER_FFTPREP: advance preparation for smoothing by FFT
% Syntax:
%   fftinfo = register_fftprep(im_sz,sigma)
% where
%   im_sz is the size of the fixed image;
%   sigma is the vector of standard deviations used for gaussian
%     smoothing.
%
% On output, fftinfo is a structure containing several fields that can be
% used by REGISTER_SMGRAD.
%
% See also: REGISTER_SMGRAD.
  
% Copyright 2006 by Timothy E. Holy
  
  % Calculate padded size (so that we don't get wrap around effects)
  tot_sz = im_sz + 3*sigma;
  fft_sz = 2.^ceil(log2(tot_sz)); % if this is huge, beware!
  fftinfo.fft_sz = fft_sz;
  h = gaussian_filter(fft_sz,sigma,true);
  fftinfo.hfft = fftn(h);
  % Set up reflection padding
  first_sz = ceil((fft_sz-im_sz)/2);
  last_sz = fft_sz-im_sz-first_sz;
  n_dims = length(im_sz);
  x = cell(1,n_dims);
  for dimIndex = 1:n_dims
    x{dimIndex} = [first_sz(dimIndex):-1:1,1:im_sz(dimIndex), ...
		   im_sz(dimIndex):-1:im_sz(dimIndex)-last_sz(dimIndex)+1];
       x{dimIndex}(x{dimIndex}<1) = 1;
       x{dimIndex}(x{dimIndex}>im_sz(dimIndex)) = im_sz(dimIndex);
  end
  fftinfo.padcoords = x;  % this is how you insert the data for reflection bdry conditions
  % Set up coords for extracting data out again (you can also use this for
  % zero-padded boundaries)
  for dimIndex = 1:n_dims
    x{dimIndex} = (1:im_sz(dimIndex))+first_sz(dimIndex);
  end
  fftinfo.snipcoords = x;  % this is how you snip out the smoothed data again
