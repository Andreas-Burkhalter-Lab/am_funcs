function data = mg_prepare_data(data,schedule)
% MG_PREPARE_DATA: generate coarse versions of data for multigrid
%
% This is a convenience function, useful especially when optimizations
% depend upon several items of data. Otherwise, it may be easier to do
% this work manually as you can inspect the size and amount of detail in
% the data.
%
% Syntax:
%   mgdata = mg_prepare_data(data,schedule)
% where
%   data: an array or a 1-by-n_data cell array of arrays
%   schedule: an (n_levels-1)-by-n_dims logical array, indicating the
%     dimensions to restrict at each transition down the sequence of
%     grids (i.e., true indicates that the given dimension should be
%     coarsened when going to the next grid).
%     Note if any of the data arrays have more than n_dims, those
%     dimensions are not subject to restriction (i.e. the schedule will
%     be set to false for those dimensions).
% and
%   mgdata is an n_levels-by-n_data cell array of arrays, with
%     mgdata{1,1} being the fine-scale version of the first data array
%     and mgdata{n_levels,1} being the coarsest version of that data
%     array.
%
% mgdata is in a format suitable for MGOPT.
%
% See also: MGOPT.
  
% Copyright 2009 by Timothy E. Holy

  if ~iscell(data)
    data = {data};
  end
  n_data = length(data);
  [n_levels,n_dims] = size(schedule);
  n_levels = n_levels+1;
  
  % Compute the coarse-resolution data (needs to be done only once) and
  % the size of the grid at each resolution
  data{n_levels,1} = [];  % to allocate space for coarse-grid versions
  for lvlIndex = 2:n_levels
    flag = logical(schedule(lvlIndex-1,:));
    for dataIndex = 1:n_data
      thisData = data{lvlIndex-1,dataIndex};
      nd = ndims(thisData);
      flagtmp = [flag repmat(false,1,nd-n_dims)];
      data{lvlIndex,dataIndex} = array_restrict(thidData,flagtmp);
    end
  end
end