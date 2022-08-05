function w = synthesize_waveform(tspike,template,trange,options)
% SYNTHESIZE_WAVEFORM: create a simulated spike recording
% Syntax:
%   w = synthesize_waveform(tspike,template,trange,options)
% where
%   tspike is a cell array, with each element containing the vector of
%     spike times (in scans) for a given "neuron" (these times are
%     unit-offset, i.e., 1 would correspond to the first scan)
%   template is a cell array, with each element a n_channels-by-
%     templatelength matrix, where
%        template{cellIndex}(channelIndex,:)
%     is the voltage waveform on the given electrode over the snippet trange
%   trange is a 2-vector, [snipstart snipend] relative to the peak of the
%     waveform (i.e., [-5 20] would mean 5 samples before the peak and 20
%     after the peak).
%   options may have the following fields:
%     noise_amplitude (default 1): the amplitude of the gaussian noise to
%       put on each channel
%     filterb, filtera: optional filtering to apply to the noise (the
%       default is to do no filtering). filtfilt is used.
% and
%   w is a n_channels-by-nscans matrix, where the total number of scans
%     is determined from the spike times and trange. The times referenced
%     in w agree with those in tspike, i.e., w(:,1) is the first scan in
%     agreement with whatever notion is in tspike, regardless of the
%     setting for trange. In particular, spike waveforms might be clipped
%     at the left edge of w.
%
% See also: SPIKESIM.
  
% Copyright 2007 by Timothy E. Holy
  
  if ~iscell(tspike)
    tspike = {tspike};
  end
  if ~iscell(template)
    template = {template};
  end
  n_cells = length(tspike);
  if (length(template) ~= n_cells)
    error(['Number of cells is not consistent between tspke and' ...
	   ' template']);
  end
  [n_channels template_length] = size(template{1});
  if (diff(trange)+1 ~= template_length)
    error('Template size does not match trange');
  end
  if (nargin < 4)
    options = struct;
  end
  options = default(options,'noise_amplitude',1);
  
  tspike_all = cat(2,tspike{:});
  trec = [0 max(tspike_all)] + trange;
  nscans = diff(trec)+1;
  spikelen = diff(trange); % off by 1, but we'll use 0 indexing below
  
  % Start with noise
  w = options.noise_amplitude * randn(n_channels,nscans);
  if isfield(options,'filterb')
    w = filtfilt(filterb,filtera,w')';
  end
  
  % Add in spikes
  for cellIndex = 1:n_cells
    thistime = tspike{cellIndex};
    thistemplate = template{cellIndex};
    for spikeIndex = 1:length(thistime)
      trng = thistime(spikeIndex):thistime(spikeIndex)+spikelen;
      w(:,trng) = w(:,trng) + thistemplate;
    end
  end
  
  % Clip the waveform to positive times
  w = w(:,max(-trange(1),0)+1:end);
  