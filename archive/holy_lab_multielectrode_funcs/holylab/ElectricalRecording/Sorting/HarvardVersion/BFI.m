function [filters,subrange,sv,wave] = BFI(spikefiles,noisefiles,nspikes,nnoise,channel,options)
% BFI: Build Filters Interactively
subset = 0;
if (nargin > 5)
  if isfield(options,'snipindx')
    subset = 1;
    snipindx = options.snipindx;
  end
else
  options = struct;
end
% First figure out how many snippets we have/channel in each file
if (subset)
  [chans,mnsnips,ssniprange] = GetSnipNums(spikefiles,options);  %This is here solely to set ssniprange
  for i = 1:length(snipindx)
    nsnips(i) = length(snipindx{i});
  end
else
  [chans,mnsnips,ssniprange] = GetSnipNums(spikefiles,options);
  chindx = find(chans == channel);
  nsnips = mnsnips(chindx,:);
end
% Now load a representative sample of the snippets for filter construction
totsnips = sum(nsnips);
fracspike = min(nspikes/totsnips,1);
rangespike = BuildRangeMF(nsnips,fracspike);
if (subset)
  indexspike = BuildIndexMF(rangespike,snipindx);
else
  indexspike = BuildIndexMF(rangespike);
end
spikes = LoadIndexSnippetsMF(spikefiles,channel,indexspike,options);
% Do the same for noise snippets. No indexing necessary here!
if ~isempty(noisefiles)
  [chans,mnoise,nsniprange] = GetSnipNums(noisefiles,options);
  chindx = find(chans == channel);
  nnoise = mnoise(chindx,:);
  totnoise = sum(nnoise);
  fracnoise = min(nnoise/totnoise,1);
  if (fracnoise*nnoise < ssniprange(2)-ssniprange(1))
    error('Do not have enough noise snippets on this channel to build filters!');
  end
  indexnoise = BuildIndexMF(BuildRangeMF(nnoise,fracnoise));
  noise = LoadIndexSnippetsMF(noisefiles,channel,indexnoise,options);
else
  noise = zeros(size(spikes,1),0);
end
alldone = 0;
hfig = ChooseWaveforms(spikes,ssniprange);
while (alldone == 0)
  waitfor(hfig,'UserData','done');
  if (~ishandle(hfig))
    warning('Operation cancelled by user');
    filters = [];
    subrange = [];
    sv = [];
    wave = [];
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
%noise = 50*randn(ssniprange(2)-ssniprange(1)+1,nnoise);
subrange = [sniprange(1)-ssniprange(1)+1, sniprange(2)-ssniprange(1)+1];
[filters,wave,sv] = Build2Filters(spikes(subrange(1):subrange(2),goodspikes),noise(subrange(1):subrange(2),:));
%[spikesfilt,IKeepFilt] = AlignSpikesFilt(spikes(:,goodspikes),filters(:,1));
%figure
%plot(spikesfilt)
%title('Aligned')
%filters = Build2Filters(spikesfilt,noise(subrange(1):subrange(2),:));
