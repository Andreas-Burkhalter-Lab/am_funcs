function [twhiso,twhisIndexo,snips,outdata] = whistimes(sng,header,options)
% WHISTIMES: determine the time intervals containing whistles
%
% This function detects whistles based on their "spectral purity" and
% mean frequency.  In a given time slice, spectral purity is defined as
% the ratio of powers between the frequency bin with the maximum power,
% and the sum of all power in that time slice.  The mean frequency is
% weighted by the power as a function of frequency.
%
% Syntax:
%   twhis = whistimes(filename,options)   (file must be a sparse sng file)
%   twhis = whistimes(sng,header,options) (to provide a sparse sng from memory)
%   [twhis,data] = whistimes(...)
%   whistimes(...)                        (to simply plot the results)
% where
%   filename is a string containing the name of a sparse sonogram file
%   sng is a sparse matrix sonogram
%   header is the header structure read by READAIHEADER
%   options is structure which may contain the following fields
%     (See also WHISTIMESDEFAULTS):
%     log: if true, computations are done with the log-power in
%       each time/frequency bin, rather than the power in that bin.
%     filterduration: the duration of the median filter (in seconds) used
%       to smooth the spectral purity and mean frequency.  If absent, no
%       filtering is performed.
%     puritythresh: the minimum spectral purity for a candidate whistle.
%       Must be between 0 and 1.
%     specdiscthresh (optional): the maximum spectral discontinuity allowed.
%       Basically, this quantity is bigger if the frequency profile changes
%       from time slice to time slice. See SPECDISCONT.  If absent, this
%       criterion is not used.
%     meanfreqthresh (optional): the minimum mean frequency (in Hz) allowed for a
%       whistle.  Noise sources have most of their power at lower
%       frequencies, even if they "bleed" into the band occupied by
%       whistles. If absent, this criterion is not used.
%     durationthresh: the minimum allowed duration of a whistle (in
%       seconds).  Events which satisfy the purity and mean frequency
%       thresholds must do so for at least a time specified by this
%       parameter.
%     mergeclose: after whistle detection is performed, merge together
%       whistles if they are separated by only a small gap, of duration
%       shorter than this parameter (in seconds).
%
% and
%   twhis is a 2-by-nwhistles matrix of starting and ending times for the
%     whistles (in seconds).
%   data is a structure array of name, value pairs; these pairs are the
%     different statistics used to assess the presence of a whistle (e.g.,
%     purity, meanfreq, specdisc).
%
% See also: WHISTIMESDEFAULTS, WHISTIMESPLOT, SPSNGPLOT, WHISSHOWPLAY, SPECPURITY, SPECDISCONT.
  
% Copyright 2002 by Timothy E. Holy <holy@pcg.wustl.edu>
% Changelog:
%    2004-08-11: Added spectral discontinuity test; make meanfreq
%    optional
%    2004-08-29: Split plotting into a separate function

% Argument parsing
  if (nargin == 2)
    options = header;
    filename = sng;
    [sng,header] = ReadSonogram(filename);
  end
  t = linspace(0,header.nscans/header.scanrate,header.columnTotal);
  % Convert physical units into discrete units
  dt = t(end)/header.columnTotal;
  filterduration = 1;
  if isfield(options,'filterduration')
    filterduration = round(options.filterduration/dt);
  end
  durationthresh = round(options.durationthresh/dt);
  df = header.scanrate/2/(header.nfreq-1); % 8/7 corrected this. Should be samprate/nfft
  % You can confirm that I'm (Mike) right because if you do:
  % f = linspace(0,header.scanrate/2,header.nfreq);
  % and look at f(end) - f(end-1) this ==
  % header.scanrate/(2*(header.nfreq-1))
  
   mag = abs(sng);
  [i,j,s] = find(mag);
  pow = mag.^2;
  if (isfield(options,'log') && options.log)
    pow = sparse(i,j,log(s/header.threshold),size(sng,1),size(sng,2));
  end
  totPower = sum(pow);
  [maxPower,~] = max(pow);
  nzindex = find(totPower); % Where all the non-zero power entries are
  specpurity = zeros(size(t));
  specpurity(nzindex) = maxPower(nzindex)./totPower(nzindex); 
  %rather than compute for every column, he just looks at the nonzero
  if isfield(options,'meanfreqthresh')
    %weight average per column so that max power has most effect on
    %meanfreq calc
    spfreq = sparse(i,j,i,size(sng,1),size(sng,2));
    meanfreq = sum(spfreq.*pow)*df;
    meanfreq(nzindex) = meanfreq(nzindex)./totPower(nzindex);
  end
  if isfield(options,'specdiscthresh')
    sd = specdiscont(pow);
  end
  if (isfield(options,'filterduration'))
    specpurity = medfilt1(specpurity,filterduration);
    if isfield(options,'meanfreqthresh')
      meanfreq = medfilt1(meanfreq,filterduration);
    end
    if isfield(options,'specdiscthresh')
      sd = medfilt1(sd,filterduration);
    end
  end
  
  % Pick out the whistles
  badindx = specpurity < options.puritythresh;
  if isfield(options,'meanfreqthresh')
    badindx(end+1,:) = meanfreq < options.meanfreqthresh;
  end
  if isfield(options,'specdiscthresh')
    badindx(end+1,:) = sd > options.specdiscthresh;
  end
  
  badindx = find(max(badindx));
  longindx = find(diff(badindx) > durationthresh);
  twhis = [t(badindx(longindx)+1); t(badindx(longindx+1)-1)];
  twhisIndex = [badindx(longindx)+1;badindx(longindx+1)-1];
  % Mike: 6/6/2014: Added +1 and -1. Really it's the zeros in between the
  % badindx values right? This is why when I calculate peak frequency
  % sometimes the ends are weird - they shouldn't really be include in the
  % snips.

  % Merge pairs separated by only a small gap
  if isfield(options,'mergeclose')
    dt = twhis(1,2:end)-twhis(2,1:end-1);
    closei = find(dt <  options.mergeclose);
    for i = length(closei):-1:1
      twhis(2,closei(i)) = twhis(2,closei(i)+1);
      twhis(:,closei(i)+1) = [];
      twhisIndex(2,closei(i)) = twhisIndex(2,closei(i)+1);
      twhisIndex(:,closei(i)+1) = [];
    end
  end
  
  % Output values, or plot?
  featurevecs(1).name = 'purity';
  featurevecs(1).value = specpurity;
  featurevecs(1).thresh = options.puritythresh;
  if isfield(options,'meanfreqthresh')
    featurevecs(end+1).name = 'meanfreq';
    featurevecs(end).value = meanfreq;
    featurevecs(end).thresh = options.meanfreqthresh;
  end
  if isfield(options,'specdiscthresh')
    featurevecs(end+1).name = 'specdisc';
    featurevecs(end).value = sd;
    featurevecs(end).thresh = options.specdiscthresh;
  end
  if (nargout > 0)
    twhiso = twhis;
    twhisIndexo = twhisIndex;
    %if nargout > 1
    %  outdata = featurevecs;
    %end
    
    % Snips!
    snips = cell(1,size(twhisIndexo,2));
    for i = 1:size(twhisIndexo,2)
        snips{i} = sng(:,twhisIndexo(1,i):twhisIndex(2,i));
    end
  else
    whistimesplot(sng,header,twhis,featurevecs);
  end
