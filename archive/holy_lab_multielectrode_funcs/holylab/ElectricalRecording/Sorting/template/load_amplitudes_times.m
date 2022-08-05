function [amplitudes_all,times_all,residual_keys, headers] = load_amplitudes_times(fitfiles,mode)
% LOAD_AMPLITUDES_TIMES: load spike times and component amplitudes
% Syntax:
%   [amplitudes,times,residual_keys,headers] = load_amplitudes_times(fitfiles,'amplitudes')
%   [minmax,times,residual_keys,headers] = load_amplitudes_times(fitfiles,'minmax')
% where
%   fitfiles is a cell array of ".fitcomp" files (produced by
%     fit_components)
% and
%   amplitudes is a cell array, each element containing the d-by-N matrix
%     of component amplitudes. (N = # of spikes in the file)
%   minmax is cell array, each containing an nchan-by-2-by-N array of
%     min/max pairs for each (de-overlapped) snippet
%   times is a cell array, each element a 1-by-N vector containing the
%     spike time (in scans)
%   residual_keys is the list of identifying keys for the residual .merec
%     file (used to insure that the components and residuals match).
%
% See also: FIT_COMPONENTS.
  
% Copyright 2008 by Timothy E. Holy
  
  if ischar(fitfiles)
    fitfiles = {fitfiles};
  end
  n_files = length(fitfiles);

  % Load the fitting amplitudes and spike times
  amplitudes_all = cell(1,n_files);
  times_all = cell(1,n_files);
  residual_keys = nan(1,n_files);
  headers = cell(1,n_files);
  for fileIndex = 1:n_files
    s = load('-mat',fitfiles{fileIndex});
    if isfield(s,'residual_key')
      residual_keys(fileIndex) = s.residual_key;
    end
    switch mode
      case 'amplitudes'
       amplitudes_all{fileIndex} = s.amplitudes;
     case 'minmax'
      amplitudes_all{fileIndex} = s.minmax;
     otherwise
      error('mode not recognized');
    end
    % Read the spike times
    times_all{fileIndex} = s.spiketimes;
    % Check that the file exists, and if not try the parent directory
    if ~exist(s.fileToFit,'file')
        [pth,basename,ext] = fileparts(s.fileToFit);
        fileTry = ['../' basename ext];
        if exist(fileTry,'file')
            warning('sorting:fileChanged','Full path %s replaced with %s',s.fileToFit,fileTry);
            s.fileToFit = fileTry;
        end
    end
    % Get the header
    headers{fileIndex} = readheader(s.fileToFit);
  end
end
