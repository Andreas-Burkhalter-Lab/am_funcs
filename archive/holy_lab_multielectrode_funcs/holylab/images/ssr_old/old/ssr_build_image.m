function [im,coverage] = ssr_build_image(features,stackIndex,Stransform)
  n_features = length(features.S);
  stack_dims = length(features.stacksize);
  if (stack_dims == 1)
    % Handle weirdness about 1d
    stacksize = [1 features.stacksize];
  else
    stacksize = features.stacksize;
  end
  im = zeros(stacksize);
  doCoverage = false;
  if (nargout > 1)
    doCoverage = true;
  end
  if doCoverage
    coverage = zeros(size(im));
  end
  coords_I = cell(1,stack_dims);
  for featureIndex = 1:n_features
    for k = 1:stack_dims
      coords_I{k} = Stransform.xmin(k,featureIndex):...
        Stransform.xmax(k,featureIndex);
    end
    im(coords_I{:}) = im(coords_I{:}) + ...
      Stransform.imT{featureIndex}*features.T(featureIndex,stackIndex);
    if doCoverage
      coverage(coords_I{:}) = coverage(coords_I{:}) + ...
        Stransform.covT{featureIndex};
    end
  end
