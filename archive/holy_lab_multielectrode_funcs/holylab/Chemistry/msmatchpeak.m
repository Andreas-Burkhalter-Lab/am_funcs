function compound = msmatchpeak(msdata,peak,options)
% MSMATCHPEAK: identify a peak (compound) in mass spec data across fractions
% Syntax:
%   compound = msmatchpeak(msdata,peak)
%   compound = msmatchpeak(msdata,peak,options)
% where
%   msdata is a structure array containing the mass spec data (see
%     MSLOAD);
%   peak is a peak structure array in one fraction (see MS2PEAK);
%   options is an option structure with the following fields:
%     thresh (default 15) is a parameter controlling how close m/z values
%       have to be in order to consider the peak as identified (see
%       below);
%     useraw (default false): if true, re-opens the raw mass spectrum
%       files and extracts all peaks within "thresh" (see above) standard
%       deviations of the mean m/z for this peak;
% and
%   compound is a structure with the following fields:
%     label: the label assigned to this compound (mean m/z from peak);
%     fraction: a vector of fraction numbers;
%     discrep: a vector of scores indicating the degree of overlap in m/z
%       with the reference peak---there is one score for each fraction
%       (bigger means less overlap, more discrepancy);
%     meanmag: a vector, containing its average magnitude in each fraction;
%     stdmag: a vector, containing the standard deviation of the magnitude
%       for each fraction;
%     ntrials: a vector, indicating the number of ms trials used for each
%       fraction;
%     peakindx: an index vector into the peak structure of each
%       fraction, indicating the peak number corresponding to this
%       compound.
%
% Perhaps you'll want to follow-up a call to this function with actually
% labelling the peaks in the data (so they're removed from future
% consideration), by calling MSLABELPEAKS.
%
% The discrepancy score is the difference between m/z means divided by the
% pooled standard deviations.  It's like Student's t, except standard
% deviations rather than standard errors.
%
% In "useraw" mode, the compound.meanmag is the average across all peaks
% satisfying the criterion.  compound.peakindx will be set to NaN, and a
% new field "raw" will be added, a cell array holding information about
% the different fractions, each of which is a cell array corresponding to
% criterion-meeting peaks in each trial.  Also, compound.ntotpeaks will
% be a vector containing the total number of peaks within the criterion
% distance in each fraction.
%
% Note: for now, this function ignores labels; only MSCHOOSEPEAK uses
%   label-avoidance.
%
% See also: MSLOAD, MS2PEAK, MSCHOOSEPEAK.
  
  if (nargin < 3)
    options = struct;
  end
  
  if ~isfield(options,'thresh')
    options.thresh = 15;
  end
  if ~isfield(options,'useraw')
    options.useraw = 0;
  end

  compound.label = mean(peak.mz);
  compound.fraction = [msdata.fraction];
  cmom2 = sum((peak.mz - compound.label).^2);
  NA = length(peak.mz);
  
  nfractions = length(msdata);
  for i = 1:nfractions
    if ~options.useraw
% $$$       npeaks = length(msdata(i).peaks);
% $$$       d = zeros(1,npeaks);
% $$$       for j = 1:npeaks
% $$$         tmp = msdata(i).peaks(j).mz;
% $$$         NB = length(tmp);
% $$$         sD = sqrt((cmom2 + sum((tmp-mean(tmp)).^2))/(NA+NB-2));
% $$$         d(j) = abs(mean(tmp) - compound.label)/sD;
% $$$       end
% $$$       [compound.discrep(i),pitmp] = min(d);
      % Find closest peak
      d = abs([msdata(i).peaks.meanmz] - compound.label);
      [dmin,pitmp] = min(d);
      tmp = msdata(i).peaks(pitmp).mz;
      NB = length(tmp);
      sD = sqrt((cmom2 + sum((tmp-mean(tmp)).^2))/(NA+NB-2));
      compound.discrep(i) = abs(mean(tmp) - compound.label)/sD;
      if (compound.discrep(i) < options.thresh)
        compound.peakindx(i) = pitmp;
        magtmp = msdata(i).peaks(pitmp).mag;
        compound.meanmag(i) = mean(magtmp);
        compound.stdmag(i) = std(magtmp);
        compound.ntrials(i) = length(magtmp);
      else
        compound.peakindx(i) = NaN;
        compound.meanmag(i) = 0;
        compound.stdmag(i) = 0;
        compound.ntrials(i) = 0;
      end
    else  % useraw: reload raw spectrum data
      m = dlmread(msdata(i).filename);
      % Split m into trials
      dm = diff(m(:,2));
      ibreak = [1;find(dm < 0)+1;length(dm)+2];
      ntrials = length(ibreak)-1;
      compound.raw{i} = cell(1,ntrials);
      for j = 1:ntrials
        mtemp = m(ibreak(j):ibreak(j+1)-1,:);
        % Identify peaks within criterion range
        peakkeep = find(abs(mtemp(:,2) - compound.label) < ...
                     sqrt(cvar)*options.thresh);
        compound.raw{i}{j} = mtemp(peakkeep,:);
      end
      allpeaks = cat(1,compound.raw{i}{:});
      compound.peakindx(i) = NaN;
      compound.meanmag(i) = mean(allpeaks(:,1));
      compound.stdmag(i) = std(allpeaks(:,1));
      compound.ntrials(i) = ntrials;
      compound.ntotpeaks(i) = size(allpeaks,1);
    end
  end