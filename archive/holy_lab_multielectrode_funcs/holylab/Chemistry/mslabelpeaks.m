function msdata = mslabelpeaks(msdata,compound)
% MSLABELPEAKS: apply compound labels to peaks
% Syntax:
%   msdata = mslabelpeaks(msdata,compound)
% where
%   msdata is a mass spec data structure
%   and compound is a compound structure such as that returned by
%     MSMATCHPEAK.
%
% See also: MSMATCHPEAK.
  
  for i = 1:length(msdata);
    if ~isnan(compound.peakindx(i))
      msdata(i).peaks(compound.peakindx(i)).label = compound.label;
    end
  end
  