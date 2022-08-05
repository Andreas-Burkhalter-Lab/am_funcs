function components_from_merec(merecFilename,chans,filenameOut,options)
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
memm = merecmm(merecFilename,'contiguous',true);
waveformRange = max(1,options.waveformRange(1));
waveformRange(2) = min(options.waveformRange(2),memm.nscans);
waveformRange(2) = 2^floor(log2(waveformRange(2))); % insure it's a power of 2
fprintf('Loading data...');
v = memm(chans(:),waveformRange);
fprintf('done\n');

fprintf('Analyzing noise parameters...');
[d,N] = size(v);
medsz = min(N,100000);
medv = median(v(:,1:medsz),2);
v = v - repmat(medv,1,N);
vnoise = mean(abs(v),2); % robust estimation of noise

thresh = 6*vnoise;
fprintf('done\n');

options = default(options,'thresh',thresh,'skip_baselineshift',true,'n_components',ceil(1.3*numel(chans)));

% If user hasn't picked a sniprange, present them with graphical data
if (~isfield(options,'sniprange') || ~isfield(options,'polarity'))
  vtmp = v(:,1:memm.scanrate);  % 1s worth of data
  [t,peakVal,peakChan] = findpeaks_multichan(vtmp,thresh,struct('polarity',0));
  sr = -200:200;
  sniplen = sr(end)-sr(1)+1;
  keepFlag = t+sr(1) > 0 & t+sr(2) < N;
  t = t(keepFlag);
  nsnips = length(t);
  snip = zeros(sniplen,nsnips);
  for i = 1:nsnips
    snip(:,i) = v(peakChan(i),t(i)+sr)';
  end
  figure
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
end

[components,lambda] = components_from_waveform(v,options);

% Massage results into expected format
s.channels = chans(:)';
s.medv = medv;
s.sniprange = options.sniprange;
s.polarity = options.polarity;
componentsc = permute(components,[2 1 3]);
sz = size(componentsc);
componentsc = reshape(componentsc,[sz(1)*sz(2) sz(3)]);
componentsc = mat2cell(componentsc,size(componentsc,1),ones(1,sz(3)));
s.templates = componentsc;
s.thresh = options.thresh;
s.merecFilename = merecFilename;
s.lambda = lambda;
if useArrayField
  s.arrayField = arrayField;
end

save('-mat',filenameOut,'-struct', 's');
