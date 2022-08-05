function b = ssr_rhs(input_data,features,options)
  filemode = iscell(input_data);
  n_features = length(features.S);
  n_stacks = length(features.stacknum);
  stack_dims = length(features.stacksize);
  % The next two variables are to handle MATLAB's weirdness about 1d
  stack_dims_md = max(2,stack_dims);
  stack_dims_sd = sum(features.stacksize > 1);
  
  % Stuff for registration
  Tinput = [eye(stack_dims_md,stack_dims_md);zeros(1,stack_dims_md)]; % empty Tform
  R = makeresampler('linear','fill');
  %R = makeresampler('linear','replicate'); % this one can be replicate even if deltaI really needs to be fill
  
  % Preallocation
  if options.calcT
    b.T = nan(n_features,n_stacks);
  end
  if options.calcS
    b.S = cell(1,n_features);
    b.pixelCount = cell(1,n_features);
    for featureIndex = 1:n_features
      b.S{featureIndex} = zeros(size(features.S{featureIndex}));
      b.pixelCount{featureIndex} = b.S{featureIndex};
    end
    im1 = ones(features.stacksize,'single');
  end
  if options.calcRegistration
    b.Registration = nan(n_features*stack_dims_sd,n_stacks);
  end
  
  % Loop over stacks
  for stackIndex = 1:n_stacks
    % Load a stack
    if filemode
      load(input_data{stackIndex});
    else
      im = input_data.smm(:,:,:,features.stacknum(i));
    end
    %if (options.calcT || options.calcRegistration)
      % Calculate the forward-transformed features
      doCoverage = true;
      Stransform = ssr_feature_transform(features,stackIndex,doCoverage);
      [im_synth,coverage] = ssr_build_image(features,stackIndex,Stransform);
      deltaim = im - im_synth;
      deltaim = deltaim .* (coverage >= 0.99);
      if options.calcRegistration
        dS = ssr_feature_gradients(Stransform.imT,stack_dims_sd);
      end
    %end
    if options.calcS
      % Extract the backwards-transformed bits of the image that overlap
      % with the features
      featureSnip = cell(1,n_features);
      pixelCount = cell(1,n_features);
      for featureIndex = 1:n_features
        offset = features.registration.offset(:,featureIndex,stackIndex)';
        % Transform the image as defined by feature's offset, cropping out
        % only the part that overlaps with the current feature.
        Tinput(stack_dims_md+1,:) = -offset; % transform image in opposite dirn
        T = maketform('affine',Tinput);
        feature_size = size(features.S{featureIndex});
        Tcrop = maketform('box',feature_size,...
          [ones(1,stack_dims_md-stack_dims) features.Xmin(:,featureIndex)+1],...
          [ones(1,stack_dims_md-stack_dims) features.Xmax(:,featureIndex)+1]);
        Ttotal = maketform('composite',fliptform(Tcrop),T);
        imT = tformarray(deltaim,Ttotal,R,1:stack_dims_md,1:stack_dims_md, ...
          feature_size,[],0);
        featureSnip{featureIndex} = imT;
        pixelCount{featureIndex} = tformarray(im1,Ttotal,R,1:stack_dims_md,1:stack_dims_md, ...
          feature_size,[],0);
      end
    end
    
    coords_I = cell(1,stack_dims);
    if options.calcT
      for featureIndex = 1:n_features
        for dimIndex = 1:stack_dims
          coords_I{dimIndex} = Stransform.xmin(dimIndex,featureIndex):...
            Stransform.xmax(dimIndex,featureIndex);
        end
        tmpProductImage = deltaim(coords_I{:}) .* Stransform.imT{featureIndex};
        b.T(featureIndex,stackIndex) = sum(tmpProductImage(:));
      end
    end
    if options.calcS
      for featureIndex = 1:n_features
        b.S{featureIndex} = b.S{featureIndex} + ...
          featureSnip{featureIndex}*features.T(featureIndex,stackIndex);
        b.pixelCount{featureIndex} = b.pixelCount{featureIndex} + pixelCount{featureIndex};
      end
    end
    if options.calcRegistration
      tmpReg = zeros(stack_dims_sd,n_features);
      for featureIndex = 1:n_features
        for dimIndex = 1:stack_dims
          coords_I{dimIndex} = Stransform.xmin(dimIndex,featureIndex):...
            Stransform.xmax(dimIndex,featureIndex);
        end
        delta_im_snip = deltaim(coords_I{:});
        for dimIndex = 1:stack_dims_sd
          tmpProductImage = dS{dimIndex,featureIndex} .* ...
            delta_im_snip;
          % The negative sign is because this has to go to the other side
          % of the equation
          tmpReg(dimIndex,featureIndex) = -sum(tmpProductImage(:)) * ...
            features.T(featureIndex,stackIndex);
        end
      end
      b.Registration(:,stackIndex) = tmpReg(:);
    end
  end % over stacks
  % Correct for the pixelCount in the spatial features
  if options.calcS
    for featureIndex = 1:n_features
      nzIndex = b.pixelCount{featureIndex} > 0;
      b.S{featureIndex}(nzIndex) = n_stacks*(b.S{featureIndex}(nzIndex) ./ b.pixelCount{featureIndex}(nzIndex));
    end
  end