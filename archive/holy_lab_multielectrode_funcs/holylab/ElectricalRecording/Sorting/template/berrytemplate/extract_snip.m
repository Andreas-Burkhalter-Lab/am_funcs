function snip=extract_snip(filename, channels, idxRawCluster, sniptimes, medv, options)
% extract_snip: extract snippets for the given raw cluster
% USAGE:
%    snip=extract_snip(filename, channels, idxRawCluster, sniptimes, medv, options)
% PRE:
%    filename, channels, idxRawCluster, sniptimes, medv: 
%       all params here are what saved by cluster_events_raw()
%    options: what saved by cluster_events_raw() plus one more field
%       maxShift which is the extra amount of samples to be read 
% POST:
%    snip: all snippets for the specified raw cluster. Each column is one snippet; 
%          each snippet is saved channel-by-channel.
  if iscell(sniptimes)
   spikeTimes=round(sniptimes{idxRawCluster});
  else
    spikeTimes = round(sniptimes);
  end
   nSnippets=length(spikeTimes);
   what=sprintf('Loading {\\\\color{blue}%d} snippets in raw cluster {\\\\color{blue}%d} ...', ...
      nSnippets, idxRawCluster);
   sniprange = options.sniprange;
   sniprange = sniprange+[-1 1]*options.maxShift; % the extra samples for alignment purpose
   snipMedian=repmat(medv, 1, diff(sniprange)+1);
   memm=merecmm(filename,'tovolts',true,'contiguous',true);
   snip = nan(length(channels)*(diff(sniprange)+1),length(spikeTimes));
   for idxSnip=1:length(spikeTimes)
      [ttSnip,memm]=memm(channels, spikeTimes(idxSnip)+sniprange);
      ttSnip=ttSnip-snipMedian;  % normalize
      ttSnip=ttSnip';
      snip(:,idxSnip)=ttSnip(:); % each column is a snippet; each snippet is saved channel-by-channel
      
      fShowProgress(idxSnip, nSnippets, what);
   end
   

function fShowProgress(idxSnip, nSnippets, what)
   if(mod(idxSnip,10)==1 || idxSnip==nSnippets)
      progress_bar(struct('progress', idxSnip, 'max', nSnippets, 'what', what, 'tex_mode', 1));
   end % if, show progress
   
   