function [M,minmax] = msdecimate(scan,options)
% options fields: n_bins_mz
  if (nargin < 2)
    options = struct;
  end
  options = default(options,'n_bins_mz',200);

  [mz,scanIndex,I] = mscollect_points(scan,options);
  
  % Convert m/z to integer values (to do binning)
  if isfield(options,'max_mz')
    mx = options.max_mz;
  else
    mx = max([scan.max_mz]);
  end
  if isfield(options,'min_mz')
    mn = options.min_mz;
  else
    mn = min([scan.min_mz]);
  end
  mzI = int32(round(((options.n_bins_mz-1)/(mx-mn)) * (mz-mn))+1);
  
  % Generate the binned intensity matrix
  M = accumarray([mzI; scanIndex]',I',[options.n_bins_mz length(scan)]);
  minmax = [mn mx];
end