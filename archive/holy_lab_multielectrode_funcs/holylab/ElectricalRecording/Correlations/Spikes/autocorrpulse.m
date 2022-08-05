function [acmean,acsem,npb,tacindx,binlist,tac,ac] = autocorrpulse(t,tbins,v)
% autocorrpulse: compute auto-correlations in an atomic signal
% 
% If a series of pulses of heights v(t_i) come at a discrete series of
% times {t_i}, compute the autocorrelation of the vs as a function of time.
%
% Syntax:
%   [acmean,acsem,npb] = autocorrpulse(t,tbins,v)
%   [acindx,binlist,npb] = autocorrpulse(t,tbins)
%   [acmean,acsem,npb,acindx,binlist,tac,ac] = autocorrpulse(t,tbins,v)
% where
%   t is a vector of times;
%   tbins is the set of right bin boundaries;
%   v is the row vector of corresponding amplitudes, or a d-by-ntimes
%     matrix of corresponding amplitudes;
% and
%   acmean is the mean of the correlation in time bins (a d-by-nbins
%     matrix);
%   acsem is the standard error in the mean across time bins;
%   npb is the number of observations per bin.
%
% If v is not supplied (or if you ask for extra output), then the output
% is:
%   acindx: a 2-by-N matrix, giving the index of pairs contributing to the
%     autocorrelation
%   binlist: a cell array of length nbins, where binlist(i) is an index for
%     the pairs (in acindx) which contribute to the ith bin;
%   npb is the number of observations per bin.
% This form is useful if the observations are such that the autocorrelation
% involves more than simple multiplication.
% 
% A final pair of optional outputs are:
%   tac, the time lag between observations;
%   ac, the product of the two observations.
%
% See also: AUTOCORRSPIKE.

% Copyright 2004 Timothy E. Holy

tmax = tbins(end);
nbins = length(tbins);
[tac,tacindx] = autocorrspike(t,tmax);
if isempty(tac)
  % No observations, so exit early
  if (nargin < 3)
    acmean = zeros(2,0);
    acsem = cell(1,nbins);
    npb = zeros(1,nbins);
  else
    d = size(v,1);
    acmean = nan(d,nbins);
    acsem = nan(d,nbins);
    npb = zeros(1,nbins);
    tacindx = zeros(2,0);
    binlist = cell(1,nbins);
    ac = nan(d,0);
  end
  return
end

% Create the lists of all pairs in a particular bin. "sort" gives an
% efficient way to do that.
[npb,binlabel] = histc(tac,[0,tbins]);
npb(end) = [];  % histc gives an extra bin for the right edge, which is never used
[sbinlabel,binindx] = sort(binlabel);
jumpindx = [find(diff(sbinlabel)),length(sbinlabel)];
binlist = cell(1,length(tbins));
startindx = 1;
for i = 1:length(jumpindx)
  binlist{sbinlabel(startindx)} = binindx(startindx:jumpindx(i));
  startindx = jumpindx(i)+1;
end

if (nargin < 3)
  acmean = tacindx;
  acsem = binlist;
  return;
end

[d,ntimes] = size(v);
if (ntimes == 1)
  v = v';    % In case v is supplied as a column vector
  [d,ntimes] = size(v);
end
if (ntimes ~= length(t))
  error('The number of observations does not match the number of times');
end

acmean = nan(d,nbins);
acsem = nan(d,nbins);

%v = v - mean(v,2);
% Compute the raw correlations
ac = v(:,tacindx(1,:)).*v(:,tacindx(2,:));

% Calculate average and sem in each time bin
for i = 1:length(binlist)
  if (npb(i))
    acmean(:,i) = mean(ac(:,binlist{i}),2);
    acsem(:,i) = std(ac(:,binlist{i}),0,2)/sqrt(npb(i));
  end
end
