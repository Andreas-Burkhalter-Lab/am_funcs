function components_from_merec(merecFilename,chans,filenameOut,options)
% COMPONENTS_FROM_MEREC: calculate waveform components from saved data
% Syntax:
%   components_from_merec(merecFilename,arrayname,filenameOut,options)
%   components_from_merec(merecFilename,channels,filenameOut,options)
% where
%  merecFilename is a string containing the name of the .merec file to
%    analyze
%  arrayname is a string, 'left' or 'right' OR
%  channels is a vector of channel numbers
%  filenameOut is the name of the file to save the data to
%  options may have the following fields:
%    waveformRange (default [1 1e7]): the maximum range of scan #s to
%      analyze (smaller is faster, but bases results on less data)
%    polarity: -1, 1, or 0, depending on whether your spikes are
%      upward-going, downward-going, or both, respectively.
%    sniprange: the range from the peak to the [beginning end] of a spike.
%  If polarity & sniprange are not known, you will be given a chance to
%  supply them based on graphical information.
%  Furthermore, any options that can be passed to COMPONENTS_FROM_WAVEFORM
%  will also be used.
%
% See also: COMPONENTS_FROM_WAVEFORM.

% Copyright 2008 by Timothy E. Holy

if (nargin < 4)
  options = struct;
end
options = default(options,'waveformRange',[1 1e7]);
if ~isfield(options,'polarity') || ~isfield(options,'sniprange')
  fprintf('Please wait. I will ask you to supply polarity and/or sniprange.\n');
end
useArrayField = false;
if ischar(chans)
  arrayField = chans;
  chans = get_hda_holylab(chans);
  useArrayField = true;
end
fu = fit_utilities;
memm = merecmm(merecFilename,'contiguous',true);
waveformRange = max(1,options.waveformRange(1));
waveformRange(2) = min(options.waveformRange(2),memm.nscans);
waveformRange(2) = 2^floor(log2(waveformRange(2))); % insure it's a power of 2
fprintf('Loading data...');
v = memm(chans(:),waveformRange);
fprintf('done\n');

fprintf('Analyzing voltage offsets...');
[d,N] = size(v);
medsz = min(N,100000);
medv = median(v(:,1:medsz),2);
v = v - repmat(medv,1,N);
fprintf('done\n');

options = default(options,'skip_baselineshift',true);

% If user hasn't picked a sniprange, present them with graphical data
if (~isfield(options,'sniprange') || ~isfield(options,'polarity'))
  vnoise = mean(abs(v),2); % robust estimation of noise
  thresh = 6*vnoise;
  n_samples_for_inspection = min(10*memm.scanrate,size(v,2));
  vtmp = v(:,1:n_samples_for_inspection);  % 10s worth of data
  [t,peakVal,peakChan] = findpeaks_multichan(vtmp,thresh,struct('polarity',0));
  sr = -200:200;
  sniplen = sr(end)-sr(1)+1;
  keepFlag = t+sr(1) > 0 & t+sr(end) < N;
  t = t(keepFlag);
  nsnips = length(t);
  snip = zeros(sniplen,nsnips);
  for i = 1:nsnips
    snip(:,i) = v(peakChan(i),t(i)+sr)';
  end
  hfig = figure;
  plot(sr,snip)
  if ~isfield(options,'polarity')
    p = input('Choose the polarity (-1, 0, or 1): ');
    if (isempty(p) || ~any(p == [-1 0 1]))
      error('Invalid input for p, quitting');
    end
    options.polarity = p;
  end
  if ~isfield(options,'sniprange')
    finalsr = input('Choose the sniprange ([min max]): ');
    if (length(finalsr) ~= 2)
      error('Invalid input for the sniprange');
    end
    options.sniprange = finalsr;
  end
  close(hfig);
end

[templates,options] = components_from_waveform(v,options);

% Massage results into expected format
s.channels = chans(:)';
s.medv = medv;
s.sniprange = options.sniprange;
s.polarity = options.polarity;
s.templates0 = templates;
[s.templates,s.lambda] = fu.svd_components(templates);
s.thresh = options.thresh;
s.merecFilename = merecFilename;
s.options = options;
if useArrayField
  s.arrayField = arrayField;
end

save('-mat',filenameOut,'-struct', 's');
