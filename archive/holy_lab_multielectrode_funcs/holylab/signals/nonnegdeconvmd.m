function [imd,convergence_info,imw] = nonnegdeconvmd(im,psfi,options)
% NONNEGDECONVMD: fast multidimensional non-negative deconvolution
% Syntax:
%    imd = nonnegdeconvmd(im,psf)
% where im is the signal you want to deconvolve (perhaps an image) and psf
% is the response function (point spread function).
% If the psf is smaller in size than the image, for the purpose of
% deconvolution it is assumed that the peak of the psf is in the center.
%
% imd is the output deconvolved signal.
%
% You can control a few parameters with the following syntax:
%    imd = nonnegdeconvmd(im,psf,options)
% where options is a structure which may have the following fields:
%    itermax (default 100): the number of conjugate-gradient steps to use;
%    tol (default 1e-3): the fractional change in error used to determine
%      convergence;
%    display (default true): show convergence progress of the algorithm.
%    noiseinfo: allows you to specify information about the noise of the
%      image, to get a better initial guess using Wiener filtering. You can
%      supply this in two forms: if supplied as a scalar, it is the rms
%      noise (in units of intensity) added to each pixel. Alternatively,
%      you can supply this as a CAMERADATA structure.
%
%   [imd,convergence_info] = nonnegdeconvmd(...)
% also returns information about the convergence of the algorithm. In this
% case, options.display defaults to false.
%
% See also: NONNEGDECONV, CAMERADATA.

% Copyright 2008 by Timothy E. Holy

  if (nargin < 3)
    options = struct;
  end
  options = default(options,'itermax',100,'tol',1e-3,'display',nargout<2,'direct',2,'pad',false);

  
  %% Check the image
  % Cast image to single, if supplied as an integer type
  cast_result = false;
  if isempty(strmatch(class(im),{'single','double'}))
    cast_result = true;
    cast_type = class(im);
    im = single(im);
  end
  if (size(im,3) == 3)
    warning('nnd:color',['You might have supplied a color image; if so,' ...
		    ' do each color channel separately']);
  end
  minim = min(im(:));
  if (minim > 0)
    warning('nnd:min_nonzero',['Your image has a non-zero minimum value.' ...
		    ' This may limit the usefulness of nonnegative' ...
		    ' deconvolution. Consider subtracting the baseline' ...
		    ' before calling this function']);
  end
  if (minim < 0)
    error(['Can''t do nonnegative deconvolution when the image has negative' ...
	   ' values.  You probably want to truncate the negatives before' ...
	   ' calling this function.']);
  end
  if (max(im(:)) <= 0)
    error('Image is zero, nothing to deconvolve');
  end
  
  %% Prepare the PSF (part 1 of 2)
  psfi = cast(psfi,class(im)); % adjust class if needed
  szpsf = size(psfi);
  onedim = false;
  if (sum(szpsf > 1) < 2)
    % PSF is along one dimension only
    onedim = true;
    transformdim = find(szpsf > 1);
    options.transformdim = transformdim;
  end
  options.onedim = onedim;

  %% Zero-pad the image up to power-of-2 size
  sz = size(im);
  sz0 = sz;
  ndims = length(sz);
  checkdims = 1:ndims;
  if onedim
    checkdims = transformdim;
  end
  sz2 = sz;
  if options.pad
    sz2(checkdims) = 2.^ceil(log2(sz(checkdims)));
  end
  if ~isequal(sz2,sz)
    sz2c = mat2cell(sz2,1,ones(1,ndims));
    im(sz2c{:}) = 0;
  end
  sz = sz2;
  
  %% Finish the preparation of the PSF
  if ~isequal(szpsf(checkdims),sz(checkdims))
    if onedim
      psf = filter_wrap(squeeze(psfi),sz(transformdim));
    else
      psf = filter_wrap(psfi,sz);
    end
  else
    psf = psfi;
    if onedim
      psf = squeeze(psf);
    end
  end
  % Normalize the PSF
  psf = psf / sum(psfi(:));
  
  %% Prepare the fourier transforms
  if onedim
    imfft = fft(im,[],transformdim);
    psffft = fft(psf);
  else
    imfft = fftn(im);
    psffft = fftn(psf);
  end
  
  %% Create the initial guess by Wiener filtering
  if isfield(options,'noiseinfo')
    [w,N2] = imwienerdeconv(imfft,psffft,options.noiseinfo);
  else
    [w,N2] = imwienerdeconv(imfft,psffft);
  end
  imwfft = imfft .* w;
  imw = flexifft(imwfft,options);
  imd = imw;
  
  if ~isequal(size(psffft),size(w))
    rep = size(w);
    rep(size(psffft) == size(w)) = 1;
    psffft = repmat(psffft,rep);
  end
  pwfft = psffft .* w;
  options.N2 = N2;
  options.onedim = onedim;

  convergence_info = [];
  if (options.direct < 2)
    [imd,convergence_info] = nnditer(imd,imwfft,pwfft,options);
  end
  if options.direct
    [imd,citmp] = nnditer(imd,imfft,psffft,options);
    if isempty(convergence_info)
      convergence_info = citmp;
    else
      convergence_info(2) = citmp;
    end
  end

  %% Output results to user
  if cast_result
    imd = cast(imd,cast_type);
  end
  % Undo the padding
  if ~isequal(sz,sz0)
    coord = cell(1,ndims);
    for i = 1:ndims
      coord{i} = 1:sz0(i);
    end
    imd = imd(coord{:});
  end
end
  
function [imd,convergence_info] = nnditer(imd,imwfft,pwfft,options)
  %% Iteratively produce the nonnegative deconvolved image
  % imd is now the Wiener-filtered image, the starting point for the
  % iterative nonnegative deconvolution.  We take a conjugate-gradient
  % approach just like in NND_CG.
  % Change-of-variables approach: let imd = y^2, and then do an
  % unconstrained conjugate-gradient minimization. Basically, since one
  % would have to do conjugate-gradient to solve any of the linear
  % problems, we might as well do those iterations on the "real" problem.
  imd(imd < 0) = 0;    % Truncate negative values
  y = sqrt(imd);
  iter = 0;
  h = 0;
  err = inf;
  deltanew = inf;
  gyv = [];
  convergence_info = struct('err',[],'converged',true);
  cpwfft = conj(pwfft);
  pwfftv = pwfft(:);
  while (1)
    % Calculate the residual (r) between the supplied image and the
    % "PSF-convolved deconvolved image"
    imdfft = flexfft(imd,options);
    rfft = pwfft.*imdfft - imwfft;
    % Calculate grad_I (the gradient with respect to intensity at each pixel)
    gIfft = cpwfft.*rfft;
    gI = flexifft(gIfft,options);
    gy = 2*(y.*gI);        % grad_y (gradient with respect to root-intensity)
    % Update the conjugate direction
    gyvOld = gyv;
    gyv = gy(:);
    deltaold = deltanew;
    deltanew = gyv'*gyv;
    if ~isempty(gyvOld)
      deltamid = gyv'*gyvOld;
    else
      deltamid = 0;
    end
    beta = max(0,(deltanew-deltamid)/deltaold); % Polak-Ribiere
    %beta = deltanew/deltaold; % Fletcher-Reeves
    h = -gy + beta*h;
    dp = h .* gy; % will be the dot product between h and gy
    if (sum(dp(:)) >= 0)
      % Oops, it isn't a descent direction, go back to using the gradient
      h = -gy;
    end
    % Minimize in the direction h. The minimum for the line search can be
    % calculated analytically, no need to do a numerical search.
    % Calculate the coefficients (c0...c4) of the quartic
    ydyfft = flexfft(y.*h,options);
    ydyfft = ydyfft(:);
    dydyfft = flexfft(h.^2,options);
    dydyfft = dydyfft(:);
    Aydyfft = pwfftv .* ydyfft;
    Adydyfft = pwfftv .* dydyfft;
    crfft = conj(rfft(:));
    cAydyfft = conj(Aydyfft);
    c0 = sum(rfft(:) .* crfft(:))/2;  % 0th order coefficient of alpha
    c1 = 2*sum(crfft.*Aydyfft);  % coefficient of 1st order term
    c2 = 2*sum(cAydyfft.*Aydyfft) + sum(crfft.*Adydyfft);   % 2nd order
    c3 = 2*sum(cAydyfft.*Adydyfft);   % coef of 3rd order term
    c4 = sum(conj(Adydyfft).*Adydyfft)/2;  % coef of 4th order term
    % Find the minimum in terms of solving for the roots of the
    % cubic. It's better not to use the cubic formulas directly, as they
    % are much more sensitive to roundoff errors than numerical
    % root-finding.
    alpha = roots(real([4*c4 3*c3 2*c2 c1]));
    % Only real roots matter
    alpha = alpha(imag(alpha) == 0);
    % Choose the one that produces the lowest value
    erralpha = real(c0 + c1*alpha + c2*alpha.^2 + c3*alpha.^3 + c4*alpha.^4);
    [minerr,minIndex] = min(erralpha);
    if options.display
      fprintf('  tot_err = %g\n',minerr);
    end
    alpha = real(alpha(minIndex));
    if (alpha < 0)
      % This hasn't happened yet, but might as well check for it
      warning('nndmd:linesearch','alpha is negative, this is weird');
    end
    % Update y
    y = y + alpha*h;
    imd = y.^2;
    errold = err;
    err = minerr/numel(imd);
    if options.display
      fprintf('Iter %d, err = %g\n',iter,err);
    end
    convergence_info.err(end+1) = err;
    if (iter >= options.itermax || abs(errold - err) < options.tol*(errold+err) || err < options.N2)
      break;
    end
    iter = iter+1;
  end
  
  if (iter >= options.itermax)
    convergence_info.converged = false;
    if options.display
      warning('nndmd:convergence','Did not converge');
    end
  end
end


function If = flexfft(I,options)
  if options.onedim
    If = fft(I,[],options.transformdim);
  else
    If = fftn(I);
  end
end

function I = flexifft(If,options)
  if options.onedim
    I = ifft(If,[],options.transformdim);
  else
    I = ifftn(If);
  end
end
