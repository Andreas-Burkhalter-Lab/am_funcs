function [u,convergence_info] = nonnegdeconv(v,template,varargin)
% NONNEGDECONV: fast one-dimensional non-negative deconvolution
% Syntax:
%    u = nonnegdeconv(v,template)
% where v and template are column vectors, v is the signal you're trying to
% match and template is the typical waveform shape.
% In cases where you have multiple inputs and multiple outputs, each column
% of v is a separate output, and template is a
% n_samples-by-n_outputs-by-n_templates matrix.  Alternatively, if each
% template is really only substantial on a few output channels, you can
% supply template as a structure array of length n_templates, where each
% element has fields chanIndex and t (the latter a matrix of size
% n_samples-by-n_active_outputs).
%
%    u = nonnegdeconv(v,template,w)
% also allows you to pass the Wiener filter associated with the template.
% This is used for creating an initial guess for u that is then iteratively
% improved by NND_CG.
%
%    u = nonnegdeconv(vfft,templatefft)
% allows you to pass the fourier transform of the input and template. In
% this case, the noise level of the appropriate w will be guessed from the
% results of Wiener filtering for you. This is therefore a particularly
% recommended approach, although you have to insure that you've done all
% the proper zero-padding, etc.
%
%   u = nonnegdeconv(...,options)
% allows you to change parameters. Here are fields that effect this code:
%    'lambdaN' (default 0): setting positive values for lambda introduces a
%      penalty on sum(u), thus tending u to zero unless "it really needs"
%      to be positive to match the data. (noise suppression)
%    'lambdaS' (default 0): coefficient of a term that enforces temporal
%      sparseness
%    's' (default []): if using sparseness, this is the vector of
%      neighbor weights
%    'display': if true, you'll see the iterative progress printed to the
%      command window (set to false, or collect 2 outputs, to turn this off)
%    'tol' (default 1e-3): the fractional change in error to decide that
%      the algorithm has converged
%    'itermax' (default 200): the maximum number of conjugate-gradient
%      steps
%
% You can also supply options that WIENER understands; these will be passed
% on, if you aren't providing w.
%
%  [u,convergence_info] = nonnegdeconv(...)
% also returns information about the convergence of the algorithm. In this
% case, options.display defaults to false.
%
% See also: NONNEGDECONVMD, WIENER, NND_CG.

% Copyright 2007 by Timothy E. Holy

  options = struct;
  have_w = false;
  for i = 1:length(varargin)
    if isstruct(varargin{i})
      options = varargin{i};
    else
      have_w = true;
      w = varargin{i};
    end
  end
  n_templates = size(template,3);
  if ~have_w
    if (~isreal(template) && ~isreal(v))
      % The fourier transforms have already been supplied. Let's assume the
      % user knows what he's/she's doing and has sized everything
      % appropriately.
      % Our task is to set the noise level.
      ct = conj(template);
      t2 = template.*ct;
      t2vec = t2(:);
      t2min = min(t2vec(t2vec > 0));
      wfft = ct./(t2 + t2min);
      u = ifft(wfft.*v);
      rat = std(u(u>0))/std(u(u<0));
      while (rat < 2)
        t2min = 100*t2min;
        wfft = ct./(t2 + t2min);
        u = ifft(wfft.*v);
        rat = std(u(u>0))/std(u(u<0));
      end
      w = wfft;
    else
      w = wiener(template,options);
    end
  end
  maxlen = max(length(v),length(w));
  if (length(template) < maxlen)
    % Pad the template with zeros
    if isreal(template)
      template(maxlen,:,:) = 0;
    else
      error('If passing fft(template), you must first resize template to the length of v');
    end
  end
  if (length(w) < maxlen)
    if isreal(w)
      w = wiener_resize(w,maxlen);
    else
      error('If passing fft(w), you must first resize w to the size of v');
    end
  end
  if (length(v) < maxlen)
    if isreal(v)
      v(maxlen,:) = 0;
    else
      error('Can''t resize v if you''re passing in fft(v)');
    end
  end
  if isreal(v)
    vfft = fft(v);
  else
    vfft = v;
  end
  if isreal(template)
    templatefft = fft(template);
  else
    templatefft = template;
  end
  if isreal(w)
    wfft = fft(w);
  else
    wfft = w;
  end
  if (nargout > 1)
    options = default(options,'display',false);
  end
  % Prepare the initial guess
  u0 = zeros([maxlen,1,n_templates],class(v));
  for templateIndex = 1:n_templates
    u0(:,1,templateIndex) = ifft(sum(wfft(:,:,templateIndex) .* vfft,2));
  end
  % Prepare filter for measuring non-sparseness (in the local temporal sense)
  options = default(options,'lambdaS',0);
  if options.lambdaS
    s = zeros(size(u0),class(u0));
    n_s = length(options.s);
    s_rep = repmat(options.lambdaS*options.s(:),[1 1 n_templates]);
    s(2:n_s+1,:,:) = s_rep;
    s(end:-1:end-n_s+1,:,:) = s_rep;
    options.lSfft = fft(s);
  end
  [u,convergence_info] = nnd_cg(vfft,templatefft,u0,options);
  if (nargout < 2 & ~convergence_info.converged)
    warning('Failed to converge');
  end
end

