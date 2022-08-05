function templates_out = parse_templates(filename,n_used)
% parse_templates: extract templates from a fine_cluster file
% This is useful if you want to use neuronal waveforms for making a
% synthetic ("fake") .merec file.
% Syntax:
%   templates_out = parse_templates(filename,n_used)
% where
%   filename is a string containing the name of the .file_cluster file
%   n_used (default 20) is the maximum number of templates to extract
%     from the file.  The n_used largest templates will be used.
% and
%   templates_out is a n_channels-by-n_scans-by-n_templates array of
%     template waveforms.
%
% See also: SYNTHESIZE_WAVEFORM.
  
% Copyright 2007 by Timothy E. Holy
  
  if (nargin < 2)
    n_used = 20;
  end
  load('-mat',filename);

  n_templates = length(templates);
  n_used = min(n_used,n_templates);
  
  % Find the amplitude distribution
  template_amplitude = nan(1,n_templates);
  for i = 1:n_templates
    template_amplitude(i) = max(abs(templates{i}));
  end
  [tas,sort_index] = sort(template_amplitude);

  % Keep the big ones
  keepIndex = sort_index(end-n_used+1:end);
  templates_out = cat(2,templates{keepIndex});
  figure
  for i = 1:n_used
    plot(templates_out(:,i))
    pause
  end
  close(gcf)
  
  % Reshape them to go back to channel,time representation
  n_channels = length(channels);
  n_scans_per_channel = length(templates{1})/n_channels;
  templates_out = reshape(templates_out,[n_scans_per_channel,n_channels,n_used]);
  templates_out = permute(templates_out,[2 1 3]);
  