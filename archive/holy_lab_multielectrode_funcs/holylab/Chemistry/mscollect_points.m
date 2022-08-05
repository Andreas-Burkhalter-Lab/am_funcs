function [mz,scanIndex,I,mzLookup] = mscollect_points(scan,options)
  % options fields: mz_max, mz_min
  if (nargin < 2)
    options = struct;
  end
  restrict_mz = isfield(options,'min_mz') || isfield(options,'max_mz');
  min_mz = 0;
  max_mz = Inf;
  options = default(options,'min_mz',min_mz,'max_mz',max_mz,'mzI',false);
  
  n_scans = length(scan);
  sIc = cell(1,n_scans);
  if ~restrict_mz
    for i = 1:n_scans
      sIc{i} = repmat(int32(i),[1 length(scan(i).mz)]);
    end
    mz = [scan.mz];
    I = [scan.intensity];
    scanIndex = cat(2,sIc{:});
  else
    mzc = cell(1,n_scans);
    Ic = cell(1,n_scans);
    for i = 1:n_scans
      flag = scan(i).mz >= options.min_mz & ...
	     scan(i).mz <= options.max_mz;
      mzc{i} = scan(i).mz(flag);
      Ic{i} = scan(i).intensity(flag);
      sIc{i} = repmat(int32(i),[1 length(mzc{i})]);
    end
    mz = cat(2,mzc{:});
    I = cat(2,Ic{:});
    scanIndex = cat(2,sIc{:});
  end
  
  if options.mzI
    % Convert the mz value to an integer
    % First let's determine the resolution
    % The resolution seems roughly log-spaced, but not precisely with a
    % base at 0. So offset it to make it a bit more accurate.
    tI = [scan.totIntensity];
    [mxtI,maxScanI] = max(tI);
    mztmp = scan(maxScanI).mz;
    mzmin = min(mz);
    offset = 3/4 * mzmin;
    lmz = log10((mztmp - offset)/(mzmin - offset));
    res = median(diff(lmz));
    % Now convert the mz value
    mz = log10((mz - offset)/(mzmin - offset));
    mz = round(mz/res)+1;
    mzmax = max(mz);
    mzLookup = 10.^((0:mzmax-1)*res);
    mzLookup = mzLookup * (mzmin-offset) + offset;
  end
end
