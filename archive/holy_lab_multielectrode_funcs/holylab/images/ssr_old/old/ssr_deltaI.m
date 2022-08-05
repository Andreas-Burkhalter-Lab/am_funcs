function [deltaI,coverage] = ssr_deltaI(im,features,stackIndex)
  doCoverage = nargout > 1;
  Stransform = ssr_feature_transform(features,stackIndex,doCoverage);
  if doCoverage
    [im_synth,coverage] = ssr_build_image(features,stackIndex,Stransform);
  else
    im_synth = ssr_build_image(features,stackIndex,Stransform);
  end
  deltaI = im - im_synth;
  return

  
  % The legacy code below has been parcelled out into ssr_feature_transform
  % and ssr_build_image
  
  deltaI = im;
  if (nargout > 1)
    coverage = zeros(size(im));
  end
  n_features = length(features.S);
  stack_dims = length(features.stacksize);
  stack_dims_md = max(2,stack_dims); % to handle MATLAB's weirdness about 1d
  Tinput = [eye(stack_dims_md,stack_dims_md);zeros(1,stack_dims_md)]; % empty Tform
  R = makeresampler('linear','fill');
  coords_I = cell(1,stack_dims);
  coords_feature = cell(1,stack_dims);
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
    xmax_crop = min(xmax,size(im)-1);
    for k = 1:stack_dims
      coords_I{k} = xmin_crop(k)+1:xmax_crop(k)+1;
      coords_feature{k} = coords_I{k} - xmin(k);
    end
    % Subtract off the transformed feature from the image
    deltaI(coords_I{:}) = deltaI(coords_I{:}) - ...
      imT(coords_feature{:})*features.T(featureIndex,stackIndex);
    if (nargout > 1)
      % Keep track of how well each pixel in the image is included. This
      % analysis acknowledges the fact that edge pixels are only partially
      % accounted-for, and indeed keeps track of the degree of coverage.
      % This version is accurate, but slower
      covT = tformarray(ones(size(features.S{featureIndex})),...
        T,R,1:stack_dims_md,1:stack_dims_md,outsz,[],0);
      coverage(coords_I{:}) = coverage(coords_I{:}) + ...
        covT(coords_feature{:});
      % This alternative version is only partially accurate, in that it
      % ignores the partial coverage of the edge pixels. But it's faster!
      %coverage(coords_I{:}) = 1;
    end
  end
  