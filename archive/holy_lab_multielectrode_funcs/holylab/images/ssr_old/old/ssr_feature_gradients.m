function dS = ssr_feature_gradients(S,stack_dims)
  % Precalculate the spatial derivatives, so that we can just look them
  % up. The only tricky bit is getting the edges right. For practical
  % reasons, it's best to have the output size be the same as the input
  % size, so we use MATLAB's "gradient" rather than diff. It's slower
  % but more technically correct
  
  % Most of the annoying stuff here is simply because MATLAB doesn't handle
  % 1d gracefully
  n_features = length(S);
  % Calculate the dimensionality. This is non-trivial because something
  % that is 1-by-n is really 1-dimensional (as far as gradient is
  % concerned) rather than 2-dimensional.
  % But since not all features will be of the same size, we have to make
  % sure we're not fooled.
%   max_feature_size = size(S{1});
%   for featureIndex = 2:n_features
%     sz = size(S{featureIndex});
%     % Insure they have the same dimensionality
%     sz(end+1:length(max_feature_size)) = 1;
%     max_feature_size(end+1:length(sz)) = 1;
%     max_feature_size = max(max_feature_size,sz);
%   end
%   max_feature_size(max_feature_size == 1) = [];
%   stack_dims = length(max_feature_size);
% $$$     colons = {':'};
% $$$     colons = repmat(colons,1,stack_dims);
% $$$     dS = cell(stack_dims,n_features);
% $$$     for i = 1:n_features
% $$$       for j = 1:stack_dims
% $$$ 	dStmp = diff(S{i},1,j);
% $$$ 	colons_j = colons;
% $$$ 	colons_j{j} = size(S{i},j);
% $$$ 	dStmp(colons_j{:}) = -Stmp(colons_j{:});
% $$$ 	dS{j,i} = dStmp;
% $$$       end
% $$$     end
  if (nargin < 2)
    stack_dims = ndims(S{1});
  end
  dS = cell(stack_dims,n_features);
  for i = 1:n_features
    if ~isempty(S{i})
      sz = size(S{i});
      sz = [sz ones(1,stack_dims-length(sz))]; % Fill with trailing 1s, where needed
      dimFlag = (sz>1);
      [dS{dimFlag,i}] = gradient(S{i});
      if any(~dimFlag)
        warning('Taking gradient of image with unit dimension');
        [dS{~dimFlag,i}] = deal(zeros(size(S{i})));
      end
      if (sum(dimFlag) > 1)
        % gradient swaps 1st & 2nd coord
        dimIndex = find(dimFlag);
        dS(dimIndex,i) = dS(dimIndex([2 1 3:length(dimIndex)]),i);
      end
    else
      [dS{:,i}] = deal(nan(size(S{i})));
    end
  end
  