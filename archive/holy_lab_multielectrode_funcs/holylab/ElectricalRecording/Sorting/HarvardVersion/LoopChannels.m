function LoopChannels(outfile,spikefiles,noisefiles,fchannels,channels,params)
% LoopChannels: shape sorting of snippet waveforms
% Two calling modes:
%        LoopChannels(outfilename,spikefiles,noisefiles,fchannels,channels,params)
%                outfilename: name of output .mat file
%                spikefiles: cell array of spike snippet filenames
%                noisefiles: cell array of noise snippet filenames
%                fchannels (optional): set of channels to use in filter
%                  construction. If unspecified, uses all. Can also pass
%                  the string 'all'.
%                channels (optional): set of channels to analyze. Same
%                  behavior as above.
%                params (optional): constants to control behavior. See
%                  DCDefParams for field names.  Can also add the following
%                  fields:
%                    machfmt: control the endian status of the opened
%                      snippet files
%
% LoopChannels(outfilename)
%
%  Use this mode when you've already gotten started on sorting, and want to
%  continue where you left off. outfilename must already exist.
  
if (nargin < 6)
  params = DCDefParams([]);
end
% Find out if output file already exists.
% If it does, load it in and start appending
matindx = findstr(outfile,'.mat');
if (isempty(matindx))
  outfull = [outfile,'.mat'];
else
  outfull = outfile;
end
fnstruct = dir(outfull);
if (isempty(fnstruct))                % If file doesn't already exist
  if (nargin < 5 | strcmp(channels,'all'))
    channels = GetSnipNums(spikefiles,params);         % Do all the channels
  end
  if (nargin < 4 | strcmp(fchannels,'all'))
    fchannels = channels;
  end
  % First, build default filters
  % Step 1: load in representatives of all channels
  fprintf('Building default filters:\n  Reading sample spike snippets from all channels...\n');
  numperfile = 5000/length(spikefiles);
  spikes = [];
  opts = params;
  opts.channels = fchannels;
  for i = 1:length(spikefiles)
    [newspikes,ssniprange] = ReadFromAllChannels(spikefiles{i},numperfile,opts);
    spikes = [spikes,newspikes];
  end
  %[spikes,ssniprange] = ReadFromAllChannels(spikefiles{10},5000);
  %spikes = AlignSpikesPeak(spikes);
  %ssniprange = ssniprange+[1,-1];
  % Step 2: do the same for noise
  if (length(noisefiles) > 0)
    numperfile = 5000/length(noisefiles);
    noise = zeros(size(spikes,1),0);
    fprintf('  Reading noise snippets from all channels...\n');
    for i = 1:length(noisefiles)
      newnoise = ReadFromAllChannels(noisefiles{i},numperfile,params);
      noise = [noise,newnoise];
    end
  else
    noise = zeros(size(spikes,1),0);
  end
  %fprintf('  Reading noise snippets from all channels in first noise file...\n');
  %noise = ReadFromAllChannels(noisefiles{1},5000);
  % Step 3: let user check over spike waveforms
  % This section is just copied from BFI
  fprintf( '  Preparing to select waveforms...\n');
  alldone = 0;
  hfig = ChooseWaveforms(spikes,ssniprange);
  while (alldone == 0)
    waitfor(hfig,'UserData','done');
    if (~ishandle(hfig))
      fprintf('Operation cancelled by user\n');
      return
    end
    goodspikes = getappdata(hfig,'GoodSpikes');
    sniprange = getappdata(hfig,'NewRange');
    if (length(goodspikes) <= sniprange(2)-sniprange(1))
      errordlg('Do not have enough spikes on this channel to build filters! Select more, or cancel.','','modal');
      set(hfig,'UserData','');
    else
      alldone = 1;
    end
  end
  close(hfig);
  % Step 4: build the filters
  subrange = [sniprange(1)-ssniprange(1)+1, sniprange(2)-ssniprange(1)+1];
  [deffilters,wave,sv] = Build2Filters(spikes(subrange(1):subrange(2),goodspikes),noise(subrange(1):subrange(2),:));
  assignin('base','DefaultFilt',deffilters);
  % Step 5: graphical output
  % This is stolen from DoChanFunctions
  figure('Name','Default filters','Position',[23   541   895   180]);
  subplot(1,3,1)
  hlines = plot(sv(1:min([15 length(sv)])),'r.');
  set(hlines,'MarkerSize',10);
  ylabel('Singular values');
  set(gca,'Tag','SVAxes');
  subplot(1,3,2)
  plot(wave);
  ylabel('Waveforms');
  set(gca,'XLim',[1 size(wave,1)]);
  set(gca,'Tag','WaveAxes');
  subplot(1,3,3)
  plot(deffilters);
  ylabel('Filters');
  set(gca,'XLim',[1 size(deffilters,1)]);
  set(gca,'Tag','FiltAxes');
  % Now get going on channel looping
  chanclust = {};
  istart = 1;
else                        % Output file already exists, let's append
fprintf(sprintf('Continuing where we left off with file %s...\n',outfile));
load(outfile)
istart = length(chanclust)+1;
params.ClustNumOffset = 0;
for k = 1:length(chanclust)
  params.ClustNumOffset = params.ClustNumOffset + size(chanclust{k},1);
end
end
warning off
for k = istart:length(channels)
  hfig = DoChannel(spikefiles,noisefiles,channels(k),deffilters,subrange,params);
  hcancel = findobj(gcf,'Tag','CancelButton');
  set(hcancel,'String','Quit');
  waitfor(hfig,'UserData','done');
  if (~ishandle(hfig))
    return
  end
  % Convert the index-based storage into spike times
  tmpindx = getappdata(hfig,'clflindx');
  t = getappdata(hfig,'t');
  params = getappdata(hfig,'params');
  scanrate = getappdata(hfig,'scanrate');        % This won't change with k
  close(hfig);
  nclust = size(tmpindx,1);
  nfiles = size(tmpindx,2);
  tclust = cell(nclust-1,nfiles);
  for i = 2:nclust                % Don't save the unassigned channel
    for j = 1:nfiles
      tclust{i-1,j} = t{j}(tmpindx{i,j});
    end
  end
  chanclust{k} = tclust;
  % Make a file copy after each channel
  % (in case something goes wrong)
  save(outfile,'chanclust','scanrate','spikefiles','noisefiles','deffilters','subrange','channels');
  % Go on to the next cell #s on the next channel
  params.ClustNumOffset = params.ClustNumOffset+nclust-1;
end
return

