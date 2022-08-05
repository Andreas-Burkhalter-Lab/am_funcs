function vno_hda_presort(options)
% VNO_HDA_PRESORT: automated component of spike sorting on high-density arrays
% Syntax: call with no arguments.
% Or, optionally, provide a structure with the following optional fields:
%     polarity (default -1): -1, 1, 0 for downward, upward, or either
%       direction spikes
%     sniprange (default [-7 50]): the range in scans relative to peak of
%       the whole spike waveform
%     normalize (default = false): normalize spikes to abs(max)
%     max_components (default = Inf): the maximum number of components used
%       in fitting waveforms.
%     skip_fitting (default = false): if true doesn't fit the waveforms to
%       templates (use this if you already have .fitcomp files)
%     skip_clustering (default = false): if true, does not cluster the fit
%       waveforms
%     skip_stimulus (default = false): if true, doesn't parse the merec
%       files to extract stimulus timing
%     skip_envelopes (default = false): if true, doesn't generate .env
%       files.
%     feedback: if present, overrides the default selection for the
%       feedback channel. (If this field is set, you can also supply other
%       options to merec2vlv).

if (nargin < 1)
  options = struct;
end
options = default(options,'polarity',-1,'sniprange',[-7 50], 'normalize', false, 'max_components', Inf, 'skip_fitting',false,'skip_clustering', false, 'skip_stimulus',false,'skip_envelopes',false);

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

% Pick a file to train templates
if (n_files == 1)
  selectedTemplates = 1;
else
  [selectedTemplates,ok] = listdlg('ListString',files_to_process,'SelectionMode','single',...
    'InitialValue',ceil(n_files/2),...
    'PromptString','Pick a file to train templates:');
  if ~ok
    return
  end
end

% Select which side(s) of the array are worth analyzing
memm = merecmm(fls{selectedTemplates},'contiguous',true);
nscans = min(memm.nscans,10*memm.scanrate);  % 10s worth of data
sides = {'left','right'};
n_sides = length(sides);
hfig = figure;
for i = 1:n_sides
  chan = get_hda_holylab(sides{i});
  v = memm(chan,[1 nscans]);
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

if ~options.skip_fitting
  % Create templates
  for i = 1:n_sides
    thisSide = sides{selected(i)};
    [success,msg] = mkdir(thisSide);
    if ~success
      error(msg)
    end
    tmpfn = [thisSide filesep thisSide '.templates'];
    if ~exist(tmpfn,'file')
      components_from_merec(fls{selectedTemplates},thisSide,...
        tmpfn,...
        struct('polarity',options.polarity,'sniprange',options.sniprange));
    else
      fprintf('%s exists, skipping creation of components\n',tmpfn);
    end
  end
  
  % Do the fitting
  for i = 1:n_sides
    thisSide = sides{selected(i)};
    cd(thisSide)
    tmplFile = [thisSide '.templates'];
    fitops = options;
    fitops.saveDirectory = '.';
    fitops.isSaveRes = true;
    for j = 1:n_files
      fit_components(tmplFile,fls{j},fitops);
    end
    cd('..')
  end
end

if ~options.skip_clustering
  % Cluster the fits
  for i = 1:n_sides
    thisSide = sides{selected(i)};
    cd(thisSide)
    compfiles = cell(1,n_files);
    for j = 1:n_files
      compfiles{j} = [basefiles{j} '.fitcomp'];
    end
    cluster_fits(compfiles,'sort1.mat',options);
    cd('..')
  end
end

% Parse the stimulus
if ~options.skip_stimulus
  fprintf('Parsing the stimulus...');
  for j = 1:n_files
    outname = [basefiles{j} '.vlv'];
    if ~exist(outname,'file')
      if isfield(options,'feedback')
        merec2vlv(outname,fls{j},options);
      else
        merec2vlv(outname,fls{j});
      end
    end
  end
  fprintf('done\n');
end

% Generate envelope files
if ~options.skip_envelopes
  fprintf('Calculating envelopes...');
  for j = 1:n_files
    outname = [basefiles{j} '.env '];
    if ~exist(outname,'file')
      [tStatus, tOutput]=system(['calenv -s 100 -o ' outname ' ' fls{j}]);
    end
  end
  fprintf('done\n');
end


