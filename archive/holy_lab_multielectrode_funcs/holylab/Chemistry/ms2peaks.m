function [peaks,rawmz] = ms2peaks(m,options)
% MS2PEAKS: extract peaks from mass spectrum data
% Syntax:
%   peaks = ms2peaks(m)
%   [peaks,rawmz] = ms2peaks(m,options)
% where
%   m is the input mass spectrum data, a matrix in which the first column
%     is magnitude and the second column is m/z; m is assumed to consist
%     of one or more trials (scans across m/z).
%   options is an optional structure with the following fields:
%     minfractrials (default 0.9): the minimum fraction of trials that
%       a given peak needs to be detected on in order to be retained;
%     mzgroup (default 0.05): tolerance on grouping together m/z values;
%     trials (default 'all'): the trial numbers to keep;
% and
%   peaks is a structure array, each element containing the data about
%     the given peak:
%       mz: the measured m/z ratio on each trial;
%       meanmz: m/z averaged across trials;
%       mag: the magnitude of the peak on each trial;
%       maxmag: the largest magnitude on any trial;
%       trial: a vector of trials numbers in which this peak was seen;
%       label: a label (scalar number) for this peak; the label is set to
%         0 in preparation for future cross-fraction labelling (see
%         MSMATCHPEAKS).
%   rawmz is the raw m/z spectrum for the trial with the largest amplitude;
%
  
% Copyright 2005 by Timothy E. Holy
  
  % Parse arguments
  if (nargin < 2)
    options = struct;
  end
  if ~isfield(options,'minfractrials')
    options.minfractrials = 0.9;
  end
  if ~isfield(options,'mzgroup')
    options.mzgroup = 0.05;
  end
  
  % Find trial breaks & split m up into trials, labelling each row by its
  % trial number
  dm = diff(m(:,2));
  ibreak = [1;find(dm < 0)+1;length(dm)+2];
  ntrials = length(ibreak)-1;
  mtrial = cell(1,ntrials);
  for i = 1:ntrials
    mtrial{i} = [m(ibreak(i):ibreak(i+1)-1,:), ...
                 repmat(i,ibreak(i+1)-ibreak(i),1)];
  end
  
  % Keep only requested trials
  if (isfield(options,'trials') && isnumeric(options.trials))
    mtrial = mtrial(options.trials);
  end
  ntrials = length(mtrial);
  
  % Glue back together (perhaps with trials eliminated)
  m = cat(1,mtrial{:});
  
  % Find clusters of peaks
  [smz,indxsort] = sort(m(:,2));   % Sort by m/z ratio
  m = m(indxsort,:);
  smzshift = meanshift1d1(smz,options.mzgroup); % Aggregate by m/z
  ipkbreak = [1; find(diff(smzshift) > options.mzgroup)+1; length(smz)+1];
                      % Look for gaps in aggregated m/z
  ncandidatepeaks = length(ipkbreak)-1;

  % For each cluster, check the number of different trials present, and
  % keep only the ones that include data from a sufficient number of trials
  options.minntrials = floor(options.minfractrials*ntrials);
  for i = 1:ncandidatepeaks
    rng = ipkbreak(i):ipkbreak(i+1)-1;
    realpeak(i) = (length(unique(m(rng,3))) >= options.minntrials);
  end
  realpeak = find(realpeak);
  npeaks = length(realpeak);
  
  % Populate the peaks structure with data
  peaks(npeaks).mz = 0;  % A sneaky method for forcing allocation
  for i = 1:npeaks
    rng = ipkbreak(realpeak(i)):ipkbreak(realpeak(i)+1)-1;
    % Check for cases where one trial contributes multiple observations,
    % and pick the observation most consistent with the median magnitude
    medmag = median(m(rng,1));
    [trialsort,tsindx] = sortrows([m(rng,3) abs(m(rng,1)-medmag)]);
      % Sort back into trial order; in cases where two observations have
      % the same trial number, break the tie by sorting them in order of
      % decreasing match to the median magnitude
    rng = rng(tsindx);
    dts = diff(trialsort(:,1));
    if any(dts == 0)
      rng(find(dts==0)+1) = [];
    end
    peaks(i).mz = m(rng,2);
    peaks(i).meanmz = mean(peaks(i).mz);
    peaks(i).mag = m(rng,1);
    peaks(i).maxmag = max(peaks(i).mag);
    peaks(i).trial = m(rng,3);
    peaks(i).label = 0;
  end
  
  % Find the biggest peak, and then identify the trial which produced the
  % biggest response to that peak. Use that trial for the "rawmz" spectrum
  mnmag = zeros(1,npeaks);
  for j = 1:npeaks
    mnmag(j) = mean(peaks(j).mag);
  end
  [mxmnr,bigpindx] = max(mnmag); % bigpindx is the index of the biggest peak
  [amp,bigtindx] = max(peaks(bigpindx).mag); % bigtindx is the trial w/
                                             % biggest amplitude
  rawmz = mtrial{bigtindx}(:,1:2);
