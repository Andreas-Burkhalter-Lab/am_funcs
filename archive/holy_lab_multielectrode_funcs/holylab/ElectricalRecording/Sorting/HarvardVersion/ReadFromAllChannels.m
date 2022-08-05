function [snipm,sniprange] = ReadFromAllChannels(filename,num,options)
fprintf('       %s\n',filename);
if (nargin > 2)
  h = readheader(filename,options);
else
  h = readheader(filename);
end
if (nargin > 2)
  [chans,chanIndex] = intersect(h.channels,options.channels);
else
  chans = h.channels;
  chanIndex = 1:length(chans);
end
totsnips = sum(h.numofsnips(chanIndex));
fracsnips = min(1,num/totsnips);
opts = options;
for i = 1:length(chans)
  opts.maxsnip = ceil(fracsnips*h.numofsnips(chanIndex(i)));
  snips{i} = LoadSnip(filename,chans(i),opts);
end
%for i = 1:length(chans)
%        if (~isempty(snips{i}))
%                figure
%                plot(snips{i});
%                title(sprintf('Channel %d',chans(i)));
%        end
%end
snipm = cat(2,snips{:});
sniprange = h.sniprange;
return
