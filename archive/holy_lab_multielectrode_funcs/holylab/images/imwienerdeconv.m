function [w,N2] = imwienerdeconv(im,psf,camera)
% IMWIENERDECONV: Wiener deconvolution of images with noise handling
% 
% This is an alternative to Matlab's DECONVWNR, with more convenient
% syntax (and probably better performance) in cases where you have a good
% understanding of the noise characteristics of your images.
%
% Syntax:
%   wfft = imwienerdeconv(im,psf)
%   wfft = imwienerdeconv(im,psf,camera)
%   wfft = imwienerdeconv(im,psf,noiseamplitude)
%   [wfft,N2] = imwienerdeconv(im,psf,camera)
% where
%   im is your image (filtered and with added noise)
%   psf is the point spread function
%   camera is a CAMERADATA structure
%   noiseamplitude is the rms noise added to each pixel of the image
%     after filtering (an alternative to the camera structure)
% and
%   wfft is the fourier transform of the Wiener filter
%   N2 is the average pixel noise (for cases where you supply the
%     camera structure)
%
% The first syntax, without any noise data, tries to infer the noise
% level by the requirement that the deconvolved image be "mostly"
% nonnegative.  This is approximate at best; in such cases you should
% also consider trying the Matlab algorithm.
  
% Copyright 2008 by Timothy E. Holy
  
  have_camera_structure = false;
  if (nargin > 2)
    if isnumeric(camera)
      sigmap = camera;
    elseif isstruct(camera)
      have_camera_structure = true;
    end
  else
    sigmap = 0;  % we don't have a noise estimate
  end
  
  datatype = class(psf);
  if ~strcmp(datatype,'single') && ~strcmp(datatype,'double')
    error('psf should be single or double');
  end
  if isreal(im)
    % Image supplied as a real image
    imfft = fftn(cast(im,datatype));
  else
    % Image supplied as the FFT
    imfft = cast(im,datatype);
  end
  n_dims = ndims(im);
  n_pix = numel(im);
  
  if (isreal(psf) && ~isequal(size(psf),size(im)))
    psf = filter_wrap(psf,size(im));
  end
  if isreal(psf)
    psffft = fftn(psf);
  else
    psffft = psf;
  end
  if (psffft(1) ~= 1)
    psffft = psffft/psffft(1); % Normalize the filter
  end
  cpsffft = conj(psffft);
  psffft2 = psffft .* cpsffft;
  
  % Handle singleton dimensions
  szpsf = size(psffft);
  average_over_dim = szpsf < size(imfft);
  averaging = any(average_over_dim);
  if any(szpsf(average_over_dim) > 1)
    error('Size mismatch between psf and image, but not in a singleton dimension');
  end
  % Prepare for replicating back up to full size
  rep = size(im);
  rep(rep == size(psffft)) = 1;
  
  if have_camera_structure
    % Estimate shot noise contribution from the average image intensity
    immean = imfft(1,1)/numel(im) - camera.bias;
    % Compute total noise power
    N2 = camera.gain*immean + camera.sigmap^2;
  else
    N2 = sigmap^2;
  end
  
  if (N2 > 0)
    % We have a real noise estimate, which we'll use to construct the
    % Wiener filter. Since we don't have access to the true
    % underlying signal, our computation of the signal-to-noise ratio has
    % to be heuristic.  Our goal is to compare the noise power to the
    % filtered signal power (which includes filtered noise power); since
    % we won't get a good estimate of the SNR in cases where the signal
    % is near the noise, we'll just kill k-values that have power near
    % the noise level.
    imfft2 = imfft .* conj(imfft);
    if averaging
      for i = length(average_over_dim):-1:1
        if average_over_dim(i)
          imfft2 = mean(imfft2,i);
        end
      end
    end
    killflag = imfft2 < 2*n_pix*N2;
    % The estimate of which frequencies to kill is noisy, since in the case
    % of zero signal we expect the power to average to the noise
    % power. Filter it over local regions in k-space, over blocks just big
    % enough that it's unlikely that all the neighbors will have happened to
    % escaped killing strictly by chance.  Choose the neighbors to be evenly
    % spaced in real physical units.
    if ~averaging
      if have_camera_structure && isfield(camera,'pixel_spacing')
        pixel_spacing = camera.pixel_spacing;
      else
        pixel_spacing = ones(1,n_dims);
      end
      l2N = log2(n_pix);
      %s = (l2N - n_dims)/sum(2*pixel_spacing);
      s = (l2N/prod(pixel_spacing))^(1/n_dims)/2;
      while (prod(2*floor(s*pixel_spacing)+1) > l2N)
        s = s / 1.2;
      end
      while (prod(2*floor(s*pixel_spacing)+1) < l2N)
        s = s*1.2;
      end
      half_width = floor(s./pixel_spacing);
      full_width = 2*half_width+1;

      h = ones(full_width,datatype);

%       fkillflag = imfilter(killflag,h);
%       killflag = fkillflag > 0;
% 
%       % Get the edges back
%       fkeepflag = imfilter(1-killflag,h);
%       keepflag = fkeepflag > 0;

      hfft = fftn(filter_wrap(h,size(killflag)));
      fkillflag = ifftn(fftn(cast(killflag,datatype)).*hfft);
      killflag = fkillflag >= 1/numel(h);
      
      % Get the edges back
      fkeepflag = ifftn(fftn(cast(1-killflag,datatype)).*hfft);
      keepflag = fkeepflag >= 1/numel(h);
    else
      % If we're averaging, we'll presume that we don't need to worry about
      % the noise.
      keepflag = ~killflag;
    end
    
    % We'll also kill all power at any frequency for which the
    % inverse filter would be infinite
    keepflag = keepflag & (psffft2 > 0);
    
    % Make it conjugate-symmetric, so we get real results
    sz = size(imfft2);
    coords_conj = cell(1,length(sz));
    for i = 1:length(sz)
      coords_conj{i} = mod(sz(i)-(1:sz(i))+1,sz(i))+1;
    end
    keepflag = keepflag | keepflag(coords_conj{:});
    
    % Compute the noise-to-signal ratio (power)
    NSR = inf(size(imfft2));
    NSR(keepflag) = N2 * psffft2(keepflag) ./ imfft2(keepflag);
    
    % Compute the Wiener filter
    w = cpsffft ./ (psffft2 + NSR);
    
    if averaging
      % Expand w back up to the size of the image
      w = repmat(w,rep);
    end
  else
    % We don't have a noise estimate. Instead, use the fact that we know
    % that the image can't be negative to get a lower bound on the noise
    NSR = psffft2(:);
    NSR = min(NSR(NSR > 0));  % Start with the minimum nonzero value
    while (1)
      w = cpsffft./(psffft2 + NSR);
      if averaging
        w = repmat(w,rep);
      end
      imd = ifftn(w .* imfft); % Wiener-deconvolved image
      imdv = imd(:);
      rat = sum(imdv(imdv > 0).^2) / sum(imdv(imdv < 0).^2);
      if (rat > 10)
        break;
      end
      %NSR = 100*NSR;
      NSR = sqrt(NSR);  % This increases towards 1 (works since PSF is normalized)
    end
  end
  