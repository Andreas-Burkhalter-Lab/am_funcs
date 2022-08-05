% History:
%   2004-08-09: 1. Increased number of real data points used when
%                  interpsnips is on (to prevent discontinuity at
%                  endpoints)
%               (RCH)

function snips = cutsnippets(wave,chanIndex,tpeak,xmax,sniprange,soptions)
  nchannels = length(chanIndex);
  if (nchannels ~= length(tpeak))
      error('Number of channels do not match');
  end
  snips = cell(1,nchannels);
  sniplength = diff(sniprange)+1;
  snipindx = sniprange(1):sniprange(2);
  snipindxInterp = (sniprange(1)-1):(sniprange(2)+1);   % if interpsnips on, will initially sample one additional point in either direction
  for i = 1:length(chanIndex)
      nsnips = length(tpeak{i});
      snips{i} = zeros(sniplength,nsnips);
      if (soptions.interpsnips & ~isempty(snips{i}))
          for j = 1:nsnips
              tempsnips = wave(chanIndex(i),snipindxInterp+tpeak{i}(j));
              tempsnips = double(tempsnips);
              snips{i}(:,j) = 0.5*xmax{i}(j)*(xmax{i}(j)-1)*tempsnips(1:end-2) - (xmax{i}(j)-1)*(xmax{i}(j)+1)*tempsnips(2:end-1) + 0.5*xmax{i}(j)*(xmax{i}(j)+1)*tempsnips(3:end);
          end
      else
          for j= 1:nsnips
              snips{i}(:,j) =  wave(chanIndex(i),snipindx+tpeak{i}(j));
          end
      end
  end
       
