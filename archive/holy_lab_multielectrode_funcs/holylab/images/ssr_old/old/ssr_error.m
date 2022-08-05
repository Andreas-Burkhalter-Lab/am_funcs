function err = ssr_error(input_data,features)
  filemode = iscell(input_data);
  n_features = length(features.S);
  n_stacks = length(features.stacknum);
  stack_dims = length(features.stacksize);
  % The next two variables are to handle MATLAB's weirdness about 1d
  stack_dims_md = max(2,stack_dims);
  stack_dims_sd = sum(features.stacksize > 1);
  
  err = 0;
  totPixels = 0;
  doCoverage = true;
  
  for stackIndex = 1:n_stacks
    % Load a stack
    if filemode
      load(input_data{stackIndex});
    else
      im = input_data.smm(:,:,:,features.stacknum(i));
    end
    Stransform = ssr_feature_transform(features,stackIndex,doCoverage);
    [im_synth,coverage] = ssr_build_image(features,stackIndex,Stransform);
    deltaim = im - im_synth;
    deltaim = deltaim.^2;
    pixelKeep = coverage(:) >= 0.999;
    err = err + sum(deltaim(pixelKeep));
    totPixels = totPixels + sum(pixelKeep);
  end
  err = err/totPixels;