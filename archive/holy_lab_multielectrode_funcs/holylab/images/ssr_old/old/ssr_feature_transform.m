function Stransform = ssr_feature_transform(features,stackIndex,doCoverage)
  if (nargin < 3)
    doCoverage = false;
  end
  n_features = length(features.S);
  stack_dims = length(features.stacksize);
  stack_dims_md = max(2,stack_dims); % to handle MATLAB's weirdness about 1d
  Tinput = [eye(stack_dims_md,stack_dims_md);zeros(1,stack_dims_md)]; % empty Tform
  R = makeresampler('linear','fill');
  coords_feature = cell(1,stack_dims);
  Stransform.xmin = nan(stack_dims,n_features);
  Stransform.xmax = nan(stack_dims,n_features);
  Stransform.imT = cell(1,n_features);
  if doCoverage
    Stransform.covT = cell(1,n_features);
  end
  
  for featureIndex = 1:n_features
    offset = features.registration.offset(:,featureIndex,stackIndex)';
    offset_floor = floor(offset);
    offset_frac = offset-offset_floor;
    % We are really only interested in the shift by fractional pixels, as
    % the integer shift can be handled using indexing
    Tinput(stack_dims_md+1,:) = offset_frac;
    T = maketform('affine',Tinput);
    outsz = size(features.S{featureIndex})+ceil(offset_frac);
    % We fill with zeros because the features are zero outside of their
    % defined region of support
    imT = tformarray(features.S{featureIndex},T,R,...
		     1:stack_dims_md,1:stack_dims_md, ...
		     outsz,[],0);
    % Determine the portion of the image that corresponds to this
    % feature/offset combination. We have to crop if it goes over the edge
    xmin = features.Xmin(:,featureIndex,:)' + offset_floor;
    xmax = features.Xmax(:,featureIndex,:)' + offset_floor + ceil(offset_frac);
    xmin_crop = max(xmin,zeros(1,stack_dims));
    xmax_crop = min(xmax,features.stacksize-1);
    for k = 1:stack_dims
      coords_feature{k} = (xmin_crop(k)+1:xmax_crop(k)+1) - xmin(k);
    end
    Stransform.xmin(:,featureIndex) = xmin_crop'+1;
    Stransform.xmax(:,featureIndex) = xmax_crop'+1;
    Stransform.imT{featureIndex} = imT(coords_feature{:});
    if doCoverage
      % Keep track of how well each pixel in the image is included. This
      % analysis acknowledges the fact that edge pixels are only partially
      % accounted-for, and indeed keeps track of the degree of coverage.
      % This approach is accurate, but slower
      covT = tformarray(ones(size(features.S{featureIndex})),...
        T,R,1:stack_dims_md,1:stack_dims_md,outsz,[],0);
      Stransform.covT{featureIndex} = covT(coords_feature{:});
    end
  end
  