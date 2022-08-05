function w = wiener(template,options)
% WIENER: Calculate one or more Wiener filters
% Syntax:
%    w = wiener(template)
% or
%    w = wiener(template,options)
% where
%   template is a column vector (for simple 1-d deconvolution), a matrix
%     (for deconvolution of multiple signals simultaneously), or a
%     n_samples-by-n_outputs-by-n_templates array (for deconvolution with
%     multiple templates).
% options may have the following fields:
%   noise2: here you can supply the mean square noise on each output channel (as
%     a vector, one entry per channel), i.e. rms.^2 or noise amplitude
%     squared.
%   SNR (default 3): signal-to-noise ratio (assumes gaussian white noise).
%     This is ignored if you supply the noise explicitly as noise2; it
%     estimates the SNR using the supplied noise2 instead.
%   wlen: the length of the Wiener filter. Default is
%     2*(SNR+1)*n_samples.
%
% If you need to use this for filtering longer data, be sure to
% look up WIENER_RESIZE, or you'll get unexpected (& undesired) behavior.
%
% See also: WIENER_RESIZE.

% Copyright 2007 by Timothy E. Holy

  if (nargin < 2)
    options = struct;
  end
  [tlen,n_outputs,n_templates] = size(template);
  templnorm = sum(template.^2,1);
  if isfield(options,'noise2')
    noise2 = options.noise2;
    if isscalar(noise2)
      noise2 = repmat(noise2,[1 n_outputs]);
    end
    % User supplied the noise; calculate the SNR, as that will be used for
    % estimating the appropriate length of the filter
    options.SNR = sqrt(max(max(templnorm,[],3) ./ (tlen*noise2)));
  end
  options = default(options,'SNR',3);
  if ~isfield(options,'noise2')
    % Define the noise on the basis of the signal power in the template(s)
    % Use the same noise on all channels, as we don't want to have outputs
    % that are simply useless (templates are small on those channels) being
    % weighted heavily.
    noise2 = mean(template.^2,1)/options.SNR^2;
    noise2 = max(noise2(:));
  end
  if isfield(options,'wlen')
    wlen = options.wlen;
  else
    wlen = 2*ceil(options.SNR+1)*tlen;
  end
  % Pad the template with zeros
  if (tlen < wlen)
    template(wlen,:,:) = 0;
  end
  % Calculate the Wiener filter
  noise2 = noise2*wlen;
  tfft = fft(template);
  if (n_outputs > 1)
    % With multiple outputs, we have to do a matrix inversion for each
    % frequency. So loop over the frequencies
    permorder = [2 3 1];
    tfft = permute(tfft,permorder);
    wfft = nan(n_outputs,n_templates,wlen);
    for k = 1:wlen
      R = diag(noise2);
      for i = 1:n_templates
        R = R + tfft(:,i,k)*tfft(:,i,k)';
      end
      wfft(:,:,k) = R\conj(tfft(:,:,k));
    end
    w = ifft(wfft,[],3);
    w = ipermute(w,permorder);
  else
    % With a single output, the "matrix" is a scalar, so we can do things
    % more efficiently
    ctfft = conj(tfft);
    denom = repmat(noise2,[wlen 1]) + sum(tfft .* ctfft,3);
    wfft = ctfft ./ repmat(denom,[1 1 n_templates]);
    w = ifft(wfft);
  end
