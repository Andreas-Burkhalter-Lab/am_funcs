function [dr,drerr] = deltar_statistics(rtrial)
% DELTAR_STATISTICS: compute mean and sem of change in firing rate
% Syntax:
%   [dr,drerr] = deltar_statistics(rtrial)
% where
%   rtrial is an nstimuli-by-ncells cell array, each element containing a
%     ntrials-by-2 matrix, with the first column the baseline firing rate
%     and the second column the response (e.g., as in the outputs
%     returned by rmax_series).
% and
%   dr and drerr are matrices of size nstimuli-by-ncells, containing the
%     mean and sem, respectively, of the firing rate change.
%
% See also: RMAX_SERIES.
  
% Copyright 2008 by Timothy E. Holy
  
  [ntags,ncells] = size(rtrial);
  dr = nan(ntags,ncells);
  drerr = nan(ntags,ncells);
  for tagI = 1:ntags
    for cellI = 1:ncells
      dr_trial = diff(rtrial{tagI,cellI},1,2);
      ntrials = length(dr_trial);
      if (ntrials > 0)
        dr(tagI,cellI) = mean(dr_trial);
        if (ntrials > 1)
          drerr(tagI,cellI) = std(dr_trial)/sqrt(ntrials);
        end
      end
    end
  end
end
  