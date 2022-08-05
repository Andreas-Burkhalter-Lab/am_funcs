function varargout = deconv_wiener_sqrt(I,psf)
% deconv_wiener_sqrt: Calculate Wiener filter for shot-noise limited images
%
% Wiener filtering assumes Gaussian noise, which is not appropriate for
% images whose noise is primarily shot noise. However, by taking a
% square-root transform, one can convert shot noise into something that is
% distributed approximately as a Gaussian.
%
% Syntax:
%   wfft = deconv_wiener_sqrt(I,psf)
% where
%   I is the image you want to filter, where intensities must be expressed
%     in quanta (e.g., electrons). I is used only (1) to determine the
%     size of the images, and (2) to determine the "typical" intensity,
%     which it does by taking the mean. If your images have a lot of black
%     regions which you don't care about, you may want to supply a "fake"
%     image that has constant intensity equal to the typical intensity you
%     care about. Or, see the multiple-output syntax below.
%   psf is the point-spread function of the imaging system. This should
%     have odd size and be centered in the middle, but does not have to be
%     equal in size to the image.
%
% On output, wfft is the Fourier transform of the Wiener filter for
% deconvolving the images. You need to implement the square-root
% transform; see the example below.
%
% Example:
% % Create the noisy filtered image
% f = fspecial('gaussian',[13 13],2);
% counts2electrons = 10;  % the "inverse gain" of the detector
% im = counts2electrons*double(imread('cameraman.tif'));  % image expressed in quanta
% imf = imfilter(im,f);  % the filtered image
% imfnoisy = poissrnd(imf);  % the noisy image
% figure; imshowsc(imfnoisy); title('Filtered, noisy image')
% % Calculate the filter
% wfft = deconv_wiener_sqrt(imfnoisy,f);
% % Deconvolve the image
% imfnoisy_tapered = edgetaper(imfnoisy,f); % optional, reduces ringing at edges
% imfn_deconv_sqrt = ifftn(fftn(sqrt(imfnoisy_tapered)).*wfft);
% % Show the deconvolved image
% figure; imshowsc(imfn_deconv_sqrt.^2); title('Deconvolved image')
%
% In this example, if you increase counts2electrons (say, to 1000) you will
% get more aggressive deconvolution, because the noise is proportionally
% smaller.
%
% An alternate syntax is the following:
%   [wnum,wdenom1,wdenom2] = deconv_wiener_sqrt(I,psf)
% In this case, I is used only for its size. You calculate wfft yourself in
% the following way:
%    wfft = wnum ./ (wdenom1 + wdenom2/Itypical)
% where Itypical (a scalar) is the typical intensity that you "care" about
% in your images, expressed in quanta. (Itypical is typical for the image,
% not the sqrt(image).)
%
% See also: deconvwnr.

% Copyright 2011 by Timothy E Holy

  if (nargout ~= 1 && nargout ~= 3)
    error('Output must have 1 or 3 arguments');
  end
  halfsize = (size(psf)-1)/2;
  if ~isequal(halfsize,round(halfsize))
    error('psf must have odd sizes in all dimensions');
  end
  if any(psf(:) < 0)
    error('psf cannot have negative entries');
  end
  index_mid = sub2ind_matrix(2*halfsize+1,halfsize);
  if (psf(index_mid) < psf(1))
    warning('filter:centered','psf should be centered in the middle, and may not be. Call fftshift?');
  end
  % Normalize the psf
  psfsum = sum(psf(:));
  if (psfsum == 0)
    error('psf cannot be zero');
  end
  psf = psf/psfsum;
  % Embed the psf into the image at full size
  n_dims = ndims(I);
  x = cell(1,n_dims);
  center_I = ceil((size(I)+1)/2);
  for i = 1:n_dims
    x{i} = center_I(i) + (-halfsize(i):halfsize(i));
  end
  psftmp = psf;
  psf = zeros(size(I),class(psftmp));
  psf(x{:}) = psftmp;
  psf = fftshift(psf);
  % Create the "deconvolved psf", a delta function at the center
  % (but then convert to corner as with fftshift)
  psfdeconv = zeros(size(I),class(psf));
  psfdeconv(1) = 1;
  % Calculate the typical intensity
  Imean = mean(I(:));
  % Calculate the shot-noise coefficient power
  N2 = numel(psf);
  % Compute the required Fourier transforms
  psf_sqrt_fft = fftn(sqrt(psf));
  psfdeconv_sqrt_fft = fftn(sqrt(psfdeconv));
  % Calculate the filter
  psf_sqrt_fft_conj = conj(psf_sqrt_fft);
  if (nargout == 1)
    varargout{1} = psfdeconv_sqrt_fft .* psf_sqrt_fft_conj ./ ...
      (psf_sqrt_fft .* psf_sqrt_fft_conj + N2/(4*Imean));
  elseif (nargout == 3)
    varargout = {psfdeconv_sqrt_fft .* psf_sqrt_fft_conj,...
      psf_sqrt_fft .* psf_sqrt_fft_conj,...
      N2/4};
  end
    
  
