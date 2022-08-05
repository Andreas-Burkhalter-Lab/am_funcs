function multichan_multiunit(options)
% MULTICHAN_MULTIUNIT: quick multiunit analysis of multichannel recordings
% Syntax: call with no arguments.
% Or, optionally, provide a structure with the following optional fields:
%     polarity (default -1): -1, 1, 0 for downward, upward, or either
%       direction spikes
%     threshold (in volts, default 0.1): the minimum amplitude of an
%       event on at least one channel 
%     skip_stimulus (default = false): if true, doesn't parse the merec
%       files to extract stimulus timing
%     skip_envelopes (default = false): if true, doesn't generate .env
%       files.
%
% See also: VNO_HDA_PRESORT.

if (nargin < 1)
  options = struct;
end
options = default(options,'polarity',-1,'skip_stimulus',false,'skip_envelopes',false,'threshold',0.1);

% Select the files
fls = UIGetFiles('*.merec','Pick files to analyze');
n_files = length(fls);
basefiles = cell(1,n_files);
files_to_process = basefiles;
for i = 1:n_files
  [pathname,basename,ext] = fileparts(fls{i});
  basefiles{i} = basename;
  files_to_process{i} = [basename ext];
end
if isempty(fls)
  return
end

selectedTemplates = ceil(n_files/2);

% Select which side(s) of the array are worth analyzing
memm = merecmm(fls{selectedTemplates},'contiguous',true);
hdr = memm.header;
nscans = min(memm.nscans,10*memm.scanrate);  % 10s worth of data
sides = {'left','right'};
n_sides = length(sides);
hfig = figure;
for i = 1:n_sides
  chan = get_hda_holylab(sides{i});
  v = memm(chan,[1 nscans]);
  medv{i} = median(v,2);
  subplot(1,n_sides,i)
  plot(v')
  title(sides{i})
end
[selected,ok] = listdlg('ListString',sides,'SelectionMode','multiple',...
  'InitialValue',1:n_sides,...
  'PromptString','Pick the sides you want to analyze');
delete(hfig)
drawnow
if ~ok
  return
end
n_sides = length(selected);
sides = sides(selected);
medv = medv(selected);

% Parse the stimulus
if ~options.skip_stimulus
  fprintf('Creating .vlv files:\n');
  for j = 1:n_files
    outname = [basefiles{j} '.vlv'];
    if ~exist(outname,'file')
      fprintf('Working on %s...\n',outname);
      merec2vlv(outname,fls{j});
    end
  end
end

% Generate envelope files
if ~options.skip_envelopes
  fprintf('Creating .env files:\n');
  for j = 1:n_files
    outname = [basefiles{j} '.env'];
    if ~exist(outname,'file')
      fprintf('Working on %s...\n',outname);
      [tStatus, tOutput]=system(['calenv -s 100 -o ' outname ' ' fls{j}]);
    end
  end
end

% Prepare for writing snippet files
n_chans = 0;
thresh = [];
channels = [];
for j = 1:n_sides
  chan_by_group{j} = get_hda_holylab(sides{j});
  channels = [channels chan_by_group{j}(:)'];
  threshtmp = medv{j};
  threshmag = options.threshold;
  if (options.polarity == 0)
    threshtmp = [threshtmp - threshmag; threshtmp + threshmag];
  else
    threshtmp = threshtmp + options.polarity*threshmag;
  end
  threshtmp = (threshtmp - hdr.scaleoff)/hdr.scalemult;
  thresh = [thresh round(threshtmp(:))'];
end
n_chans = length(channels);

default_sniphdrc = {'[snippet]',...
  'source=multichan_multiunit',...
  ['polarity=',num2str(options.polarity)],...
  ['thresh=' sprintf('%g ',thresh)],...
  sprintf('snipbeginoffset=%d',0),...
  sprintf('snipendoffset=%d',0),...
  'interptimes=0',...
  'interpsnips=0',...
  'numofsnips=',...
  'timesfpos=',...
  'snipsfpos=',...
  'finetimesfpos=',...
  'detpeaksfpos='};
default_sniphdr = sprintf('%s\n',default_sniphdrc{:});


% Do the threshold-crossing analysis
fprintf('Finding threshold-crossings:\n');
blocksize = 1e6;
for fileIndex = 1:n_files
  outname = [basefiles{fileIndex} '.ssnp'];
  if exist(outname,'file')
    continue
  end
  fprintf('Working on %s...\n',outname);
  % Prepare the header
  memm = merecmm(fls{fileIndex});
  memm.contiguous = true;
  hdr = memm.header;
  strHdr = sprintf('SNIPPET\n%s\n%s',hdr.wholeheader(7:end),default_sniphdr);
  strHdr = [strHdr sprintf('\nsnippet input file=%s\n',[basefiles{fileIndex} '.merec'])];
  strHdr = update_value(strHdr,'channel list',num2str(channels));
  [fid,msg] = fopen(outname,'w');
  if (fid < 1)
    error(['Working on file ' basefiles{fileIndex} ': ' msg]);
  end
  count = fwrite(fid,strHdr,'char');
  if (count < length(strHdr))
    error(['Error writing ' outname]);
  end
  headersize = str2double(key2value(strHdr,'header size'));
  padsize = headersize - ftell(fid);
  fwrite(fid,zeros(1,padsize),'uint8');
  if (ftell(fid) ~= headersize)
    error(['Error padding header for file ' outname]);
  end
  % Loop over blocks
  n_blocks = ceil(memm.nscans/blocksize);
  detpeaks = repmat({cell(1,n_blocks)},1,n_sides);
  peakChan = detpeaks;
  t_all = detpeaks;
  for blockIndex = 1:n_blocks
    rng = [1 blocksize]+(blockIndex-1)*blocksize;
    if (rng(2) > memm.nscans)
      rng(2) = memm.nscans;
    end
    for sideIndex = 1:n_sides
      v = memm(chan_by_group{sideIndex}(:),rng);
      [t,detpeaks{sideIndex}{blockIndex},peakChan{sideIndex}{blockIndex}] = ...
        findpeaks_multichan(v,options.threshold,options);
      t_all{sideIndex}{blockIndex} = t + rng(1);
    end
  end
  % Consolidate over blocks
  for sideIndex = 1:n_sides
    detpeaks{sideIndex} = cat(2,detpeaks{sideIndex}{:});
    peakChan{sideIndex} = chan_by_group{sideIndex}(cat(2,peakChan{sideIndex}{:}));
    t_all{sideIndex} = cat(2,t_all{sideIndex}{:});
  end
  % Consolidate over sides
  detpeaks = cat(2,detpeaks{:});
  detpeaks = round((detpeaks - hdr.scaleoff)/hdr.scalemult);
  peakChan = cat(2,peakChan{:});
  t_all = cat(2,t_all{:});
  % Write the "snippets" (really just the spike times)
  for chanIndex = 1:n_chans
    thisGroup = find(peakChan == channels(chanIndex));
    strHdr = snipfile_append_channel(fid,strHdr,channels(chanIndex), ...
      t_all(thisGroup),[],[],detpeaks(thisGroup));
  end
  fclose(fid);
  %clear(memm);
end
