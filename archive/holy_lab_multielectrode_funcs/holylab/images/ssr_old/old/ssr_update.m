function [featuresOut,err] = ssr_update(input_data,featuresIn,x,options)
  if (nargin < 4)
    options = struct;
  end
  if ~isfield(options,'safe')
    options.safe = false; % false until regularization is included
  end
  featuresOut = featuresIn;
  fn = fieldnames(x);
  if ~isempty(strmatch('T',fn,'exact'))
    % Yes, we updated the temporal component
    featuresOut.T = featuresOut.T + x.T;
  end
  if ~isempty(strmatch('S',fn,'exact'))
    % Yes, we updated the spatial component
    n_features = length(featuresOut.S);
    for featureIndex = 1:n_features
      featuresOut.S{featureIndex} = featuresOut.S{featureIndex} + x.S{featureIndex};
    end
    % We have to normalize so that the spatial component has a mean of 1
    for featureIndex = 1:n_features
      Sscale = mean(featuresOut.S{featureIndex}(:));
      featuresOut.S{featureIndex} = featuresOut.S{featureIndex} / Sscale;
      featuresOut.T(featureIndex,:) = featuresOut.T(featureIndex,:) * ...
        Sscale;
    end
  end
  if ~isempty(strmatch('offset',fn,'exact'))
    % Yes, we updated the registration
    % This one is different because it's delta(registration), not the new
    % registration. Also, the energy functional is not convex in this
    % variable, so there's no guarantee that the step will reduce the
    % error. We have to be able to control the step size to guarantee a
    % decrement.
    % Finally, in 1d there are minor issues with the dimensionality. Deal
    % with it.
    stack_dims = size(featuresIn.registration.offset,1);
    stack_dims_sd = size(x.offset,1); %1d hack
    outIndex = (1:stack_dims_sd)+(stack_dims-stack_dims_sd);
    if options.safe
      if (length(fn) == 1 && isfield(options,'errOld'))
        err0 = options.errOld;
      else
        err0 = ssr_error(input_data,featuresOut);
      end
      featuresTmp = featuresOut;
      step_size = 2;
      iter_max = 10;
      iter = 0;
      keep_trying = true;
      while (keep_trying)
        step_size = step_size/2;
        iter = iter+1;
        featuresTmp.registration.offset(outIndex,:,:) = ...
          featuresOut.registration.offset(outIndex,:,:) + ...
          step_size*x.offset;
        err = ssr_error(input_data,featuresTmp);
        keep_trying = err > err0 && iter < iter_max;
        if keep_trying
          fprintf('Registration optimization trying again, iteration %d\n',iter);
        end
      end
      if (err < err0)
        featuresOut = featuresTmp;
      else
        error('Registration update failed to converge');
      end
    else
      % We're living dangerously, don't check for convergence
      featuresOut.registration.offset(outIndex,:,:) = ...
        featuresOut.registration.offset(outIndex,:,:) + ...
        x.offset;
    end
  end
  if (nargout > 1 && ~exist('err','var'))
    err = ssr_error(input_data,featuresOut);
  end
  