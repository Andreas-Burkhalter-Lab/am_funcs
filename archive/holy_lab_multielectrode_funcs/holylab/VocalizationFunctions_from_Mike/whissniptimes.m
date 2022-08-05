function [tsnip,filtrat] = whissniptimes(ratio,t,thresh,twidth)
% WHISSNIPTIMES: determine the timespans containing vocalizations
% tsnip = whissniptimes(ratio,t,thresh,twidth)
% where
%    ratio is the ratio of max(sig)/max(ctl)
%    t is the vector of time markers for ratio
%    thresh is the threshhold for whistle detection: max(sig)/max(ctl) > thresh
%         (default: 1.5)
%    twidth is the minimum spacing between whistles, in secs (default: 0.075)
% and
%    tsnip is the 2-by-n vector of times (in seconds) of start;end for
%        the given whistle
%
% ratio is median-filtered with a filter of width twidth/4, and then
%   threshhold-crossings are sought. Once threshhold has been crossed,
%   no new snippet is activated until at least a time twidth has passed.
%
% [tsnip,filtrat] = whissniptimes(...)
%   also returns the median filtered ratio.
%
% See also: SngRatio
if (nargin < 4)
        twidth = 0.075;
end
if (nargin < 3)
        thresh = 1.5;
end
% Time bins when amplifier is saturated are indicated with a -1
% Set the value of all these to thresh
isat = find(ratio == -1);
if any(isat)
  ratio(isat) = thresh;
end

dt = t(2)-t(1);  % compute the time interval between bins
nfilt = ceil(twidth/dt/4);
blocksize = 100000;
ratiomed = medfilt1(ratio,nfilt,blocksize);
tup = find(ratiomed(1:end-1) < thresh & ratiomed(2:end) > thresh);
tdown = find(ratiomed(1:end-1) > thresh & ratiomed(2:end) < thresh);
if (~isempty(tup) & ~isempty(tdown))
        if (tdown(1) < tup(1))
                tdown = tdown(2:end);
        end
        if (tup(end) > tdown(end))
                tup = tup(1:end-1);
        end
        tsnip = [t(tup);t(tdown)];
        % Find the ones that are closer than twidth, and merge
        deltt = diff(t(tup));
        ibad = find(deltt < twidth);
        nbad = length(ibad);
        for i = nbad:-1:1
                tsnip(2,ibad(i)) = tsnip(2,ibad(i)+1);
                tsnip(:,ibad(i)+1) = [];
        end
else
  tsnip = zeros(2,0);
end
if (nargout > 1)
  filtrat = ratiomed;
end