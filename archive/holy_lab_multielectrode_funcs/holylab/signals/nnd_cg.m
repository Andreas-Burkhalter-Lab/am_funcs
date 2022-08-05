function [u,convergence_info] = nnd_cg(vfft,templatefft,u,options)
% NND_CG: the inner core of 1-d nonnegative deconvolution
% Syntax: [u,convergence_info] = nnd_cg(vfft,templatefft,u,options)
% vfft is ntimes-by-nchannels
% templatefft is ntimes-by-nchannels-by-ntemplates or a structure of length
%   ntemplates with fields chanIndex and tfft, a
%   ntimes-by-length(chanIndex) matrix. The structure format is more
%   efficient if the typical template is defined on only a few channels,
%   because one can restrict all the multiplies and add only to the
%   relevant channels.  (Note this _should_ also be better than using
%   sparse matrices, as this is such a simple form of sparse matrix.)
% u is the initial guess, ntimes-by-1-by-ntemplates
%
% See also: NONNEGDECONV, NONNEGDECONVMD.

% Copyright 2007 by Timothy E. Holy

  % Change-of-variables approach: let u = y^2, and then do an
  % unconstrained conjugate-gradient minimization. Basically, since one
  % would have to do conjugate-gradient to solve any of the linear
  % problems, we might as well do those iterations on the "real" problem.
  if (nargin < 4)
    options = struct;
  end
  options = default(options,'lambdaN',0);
  options = default(options,'itermax',200);
  options = default(options,'tol',1e-3);
  options = default(options,'display',true);
  options = default(options,'lambdaS',0);
  if isstruct(templatefft)
    n_inputs = length(templatefft);
    for templateIndex = 1:n_inputs
      templatefft(templateIndex).ctfft = conj(templatefft(templateIndex).tfft);
    end
  else
    n_inputs = size(templatefft,3);
    ctfft = conj(templatefft);
  end
  [n_samples,n_outputs] = size(vfft);
  if (size(u,3) ~= n_inputs)
    error('The initial guess for u must contain the same number of inputs as the templates');
  end
  if options.lambdaS
    if (size(options.lSfft,3) == 1 && n_inputs > 1)
      options.lSfft = repmat(options.lSfft,[1 1 n_inputs]);
    end
  end
  
  u(u < 0) = 0;    % Truncate negative values
  y = sqrt(u);
  iter = 0;
  h = 0;
  err = inf;
  deltanew = inf;
  gyv = [];
  convergence_info = struct('err',[],'dataerr',[],'converged',true);
  %figure
  %set(gca,'YLim',[-0.5 1],'NextPlot','replacechildren');
  %pause(1)
  while (1)
    u = y.^2;
    %plot(u(1:1000)); drawnow
    % Calculate the residual
    ufft = fft(u);
    if isstruct(templatefft)
      rfft = -vfft;
      for templateIndex = 1:n_templates
        cIndex = templatefft(templateIndex).chanIndex;
        rfft(:,cIndex) = rfft(:,cIndex) + repmat(ufft(:,1,templateIndex),[1 length(cIndex)]) .* ...
          templatefft(templateIndex).tfft;
      end
    else
      if (n_outputs > 1)
        ufftrep = repmat(ufft,[1 n_outputs 1]);
        rfft = sum(templatefft.*ufftrep,3) - vfft;
      else
        rfft = templatefft.*ufft - vfft;
      end
    end
    % Calculate grad_u (the gradient with respect to u)
    if isstruct(templatefft)
      gufft = zeros([n_samples 1 n_templates],class(vfft));
      for templateIndex = 1:n_templates
        gufft(:,1,templateIndex) = sum(rfft(:,templatefft(templateIndex).chanIndex) .* ...
          templatefft(templateIndex).ctfft);
      end
    else
      if (n_inputs > 1)
        rfftrep = repmat(rfft,[1 1 n_inputs]);
        gufft = sum(ctfft.*rfftrep,2);
      else
        gufft = sum(ctfft.*rfft,2);
      end
    end
    if options.lambdaS
      gufft = gufft + options.lSfft.*ufft;
    end
    gu = ifft(gufft);
    if options.lambdaN
      gu = gu+options.lambdaN;
    end
    gy = 2*(y.*gu);        % grad_y
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
    ydyfft = fft(y.*h);
    dydyfft = fft(h.^2);
    if isstruct(templatefft)
      Aydyfft = zeros([n_samples n_outputs],class(vfft));
      Adydyfft = Aydyfft;
      for templateIndex = 1:n_inputs
        chanIndex = templatefft(templateIndex).chanIndex;
        Aydyfft(:,chanIndex) = Aydyfft(:,chanIndex) + templatefft(templateIndex).tfft .* ...
          repmat(ydyfft(:,:,templateIndex),[1 length(chanIndex)]);
        Adydyfft(:,chanIndex) = Adydyfft(:,chanIndex) + templatefft(templateIndex).tfft .* ...
          repmat(dydyfft(:,:,templateIndex),[1 length(chanIndex)]);
      end
    else
      if (n_outputs > 1)
        Aydyfft = sum(templatefft .* repmat(ydyfft,[1 n_outputs 1]),3);
        Adydyfft = sum(templatefft .* repmat(dydyfft,[1 n_outputs 1]),3);
      else
        Aydyfft = templatefft .* ydyfft;
        Adydyfft = templatefft .* dydyfft;
      end
    end
    c0 = sum(sum(rfft .* conj(rfft)))/2;  % 0th order coefficient of alpha
    dataterm = c0/n_samples;   % Just the square residual, normalized for ifft
    crfft = conj(rfft);
    cAydyfft = conj(Aydyfft);
    c1 = 2*sum(sum(crfft.*Aydyfft));  % coefficient of 1st order term
    c2 = 2*sum(sum(cAydyfft.*Aydyfft)) + sum(sum(crfft.*Adydyfft));   % 2nd order
    c3 = 2*sum(sum(cAydyfft.*Adydyfft));   % coef of 3rd order term
    c4 = sum(sum(conj(Adydyfft).*Adydyfft))/2;  % coef of 4th order term
    if options.lambdaN
      % If using noise suppression, add extra terms
      % The n_samples terms are required by the normalization of the
      % fourier transform
      c0 = c0 + sum(u(:)) * n_samples;
      c1 = c1 + 2*options.lambdaN*sum(sum(y.*h))*n_samples;
      c2 = c2 + options.lambda*sum(sum(h.^2))*n_samples;
    end
    if options.lambdaS
      % If enforcing sparseness, there are yet more terms
      % Note lambdaS has already been bundled into lSfft
      % Note that doing sum(sum(...)) on each term separately saves a whole
      % bunch of adds---wait, really?
      s0 = conj(options.lSfft .* ufft);
      s1 = conj(options.lSfft .* ydyfft);
      s2 = conj(options.lSfft .* dydyfft);
      c0 = c0 + sum(sum(s0.*ufft))/2;
      c1 = c1 + sum(sum(s1.*ufft)) + sum(sum(s0.*ydyfft));
      c2 = c2 + (4*sum(sum(s1.*ydyfft)) + sum(sum(s2.*ufft)) + sum(sum(s0.*dydyfft)))/2;
      c3 = c3 + sum(sum(s2.*ydyfft)) + sum(sum(s1.*dydyfft));
      c4 = c4 + sum(sum(s2.*dydyfft))/2;
    end
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
      %fprintf('  tot_err = %g\n',minerr);
    end
    alpha = real(alpha(minIndex));
    if (alpha < 0)
      % This hasn't happened yet, but might as well check for it
      warning('nndcg:linesearch','alpha is negative, this is weird');
    end
    % Update y
    y = y + alpha*h;
    errold = err;
    err = minerr/n_samples;  % The total error
    if options.display
      fprintf('Iter %d, err = %g\n',iter,err);
    end
    convergence_info.err(end+1) = err;
    convergence_info.dataerr(end+1) = dataterm;
    if (iter > options.itermax || abs(errold - err) < options.tol*(errold+err))
      break;
    end
    iter = iter+1;
  end
  if (iter >= options.itermax)
    convergence_info.converged = false;
  end
end

