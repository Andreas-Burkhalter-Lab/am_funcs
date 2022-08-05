function [tstart] = findst_time(data,options)
% [tstart] = findst_time(data);
% this function will pool the spike times across the array to figure out
% what the offset between stimulus presentation as marked by valve opening
% and where the cells begin to respond is.  This will be independent of the
% actual cell and compound being used.
% If no ideal time is identified, the default offset time is 2.5 sec, but
% that can be changed (def_off)
%
% data: the ephys structure
% tstart: the offset time
% 
% HAA 2/11/09, modified by TEH 2011/01/10

if (nargin < 2)
  options = struct;
end
options = default(options,'binsize',0.25,'z',3,'def_off',2.4);

% pool all the tags into a big list
all_tags = cat(1,{data(:).tag});
tags = unique(all_tags);

% take all the spike times during the valve opening, combine them across
% all cells and all cycles and identify the rise time
dt = diff(data(1).scanrange)/data(1).scanrate;
toff = data(1).toffset;
for i= 1:length(tags)
  idxs = strmatch(tags(i), all_tags);
  
  tmp_t0 = cell(1,length(idxs));
  for j = 1:length(idxs) % need to rezero the times
    thisdata = data(idxs(j));
    t0 = thisdata.scanrange;
    tmp_t = cat(2,thisdata.celltimes{:}) - t0(1);
    tmp_t = tmp_t/thisdata.scanrate + thisdata.toffset;
    tmp_t0{j} = sort(tmp_t);
    dttmp = diff(thisdata.scanrange)/thisdata.scanrate;
    if (abs(dttmp-dt)/dt > 1e-3) || (thisdata.toffset ~= toff)
      error('Time ranges are not equal')
    end
  end
  
  spikes{i} = cat(2,tmp_t0{:});
end

% Bin the spike times
all_spikes = cat(2,spikes{:});
n_bins = dt/options.binsize;
trange = [0 dt] + toff;
edges = linspace(trange(1),trange(2),n_bins+1);  % this ensures all bins same size
all_hist = histc(all_spikes,edges);

% Calculate the mean and standard deviation in the pre-stimulus period
prestimflag = edges(1:end-1) <= 0;
meanr = mean(all_hist(prestimflag));
stdr = std(all_hist(prestimflag));

% Find the first post-stimulus time when the firing rate exceeds the
% threshold
post_stim_resp = all_hist(~prestimflag);
indx = find(post_stim_resp > meanr + stdr*options.z,1,'first');
if isempty(indx)
    tstart = options.def_off;
  % error('No threshold-crossing found');
end

% Back up until you get to the first time when the firing rate dips below
% the mean value
indx = find(post_stim_resp(1:indx) < meanr+stdr,1,'last');
if isempty(indx)
  indx = 1;
end

% Convert to a time
dtbin = diff(edges(1:2));
tstart = indx*dtbin;
