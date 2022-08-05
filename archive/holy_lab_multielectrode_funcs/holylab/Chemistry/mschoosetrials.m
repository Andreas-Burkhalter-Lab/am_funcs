function msdata = mschoosetrials(datain)
% MSCHOOSETRIALS: select valid mass spectra
% When mass spectra measurements are not entirely reproducible, this lets
% you choose trials in a reasonably stationary period.
% Syntax:
%   msdata = mschoosetrials(files)
%   msdata = mschoosetrials(msdatain)
% where
%   files is a cell array containing the file names of the raw mass spec
%     data (in the special case when only one file is to be loaded, you
%     may supply just the string of the name);
%   msdata/msdatain are msdata structures of the type loaded by MSLOAD.
%
% See also: MSLOAD, MSCHOOSETRIALS_GUI.
  
% Copyright 2005 by Timothy E. Holy
  msdata = [];
  if ischar(datain)
    datain = {datain};
  end
  if (iscell(datain) && ischar(datain{1}))
    % User supplied filenames
    nfiles = length(datain);
    for i = 1:nfiles
      [msdatatmp,m] = msload(datain{i});
      trials = mschoosetrials_gui(msdatatmp,m);
      if (length(trials) == 1 && trials == -1)
        return;
      end
      if ~isempty(trials)
        %msdatatmp = msct_seltrials(msdatatmp,trials);
        % Pick peaks again, using only active range of trials
        [msdatatmp.peaks,msdatatmp.rawmz] = ms2peaks(m,struct('trials',trials));
        if isempty(msdata)
          msdata = msdatatmp;
        else
          msdata(end+1) = msdatatmp;
        end
      end
    end
  elseif isstruct(datain)
    nfiles = length(datain);
    for i = 1:nfiles
      trials = mschoosetrials_gui(datain(i));
      if (length(trials) == 1 && trials == -1)
        return;
      end
      if ~isempty(trials)
        msdatatmp = msct_seltrials(datain(i),trials);
        if isempty(msdata)
          msdata = msdatatmp;
        else
          msdata(end+1) = msdatatmp;
        end
      end
    end
  else
    error('Input not recognized');
  end

function msdata = msct_seltrials(msdata,trials)
  npeaks = length(msdata.peaks);
  for i = 1:npeaks
    [ctrials,tindx] = intersect(msdata.peaks(i).trial,trials);
    msdata.peaks(i).trial = ctrials;
    msdata.peaks(i).mz = msdata.peaks(i).mz(tindx);
    msdata.peaks(i).mag = msdata.peaks(i).mag(tindx);
    msdata.peaks(i).rank = msdata.peaks(i).rank(tindx);
  end