function [val,gradout] = register_block_penalty(mismatch,uoldc,unew,barrier_mu,options)
  % register_block_penalty: combines the "data terms" (mismatch) and
  % regularization into a single objective function
  %
  % This function is called by register_block_improve, and that function is
  % a good example of how to set up the call.
  %
  % Syntax:
  %   [val,detJ,grad] = register_block_penalty(mismatch,uoldc,unew,options)
  % where
  %   mismatch is the "data penalty", a cell array where each entry
  %     describes the total square error as a function of displacement.
  %     Each element of the cell array is an array of the dimensionality of
  %     the images, where the central element corresponds to zero
  %     displacement. This is computed by register_block_mismatch and
  %     similar functions.
  %   uoldc contains the "old" deformation, arranged as a cell array with
  %     one element per spatial coordinate (i.e., a 1-by-3 cell array in 3
  %     dimensions). Each of these elements is an array of the size of the
  %     deformation grid. All zeros corresponds to the "identity" (i.e., no
  %     change) deformation.
  %   unew is the "new" deformation that will be composed with the old one.
  %     The mismatch is evaluated just with respect to the new deformation
  %     (it is to be considered in terms of an improvement on the previous
  %     warped image), but the regularization penalty works on their
  %     composition. For compatibility with optimization functions, this is
  %     represented as a pure array (often, a vector), in
  %     ugrid-then-coordinates format (in other words, when represented as
  %     an array it should be of size [gridsz n_dims]).
  %   options is a structure that is most easily set using
  %     register_phasecorr_initialize. See the example in
  %     register_block_improve for details.
  % On output,
  %   val is the value of the total penalty, combining the data terms and
  %     the regularization.
  %   detJ is the determinant of the Jacobian. This can be used to
  %     determine problematic cells (for the regularization), in terms of
  %     cells with infinite detJ.
  %   grad is the gradient of the penalty with respect to each parameter in
  %     unew.
  %
  % See also: register_block_improve, register_block_mismatch,
  % register_phasecorr_initialize.
  
  % Copyright 2011 by Timothy E. Holy
  
  
  %% Parse arguments
  gridsz = size(mismatch);
  blocksz = options.sz_spatial./gridsz;
  n_grid = numel(mismatch);
  options = default(options,'c',1);
  calc_grad = (nargout > 1);
  if isvector(unew)
    unew = reshape(unew,[gridsz options.n_dims]);
  end
  
  %% The incompressibility-Jacobianu regularization penalty
  % though compression/elongation are allowed in this code
  if options.lambda > 0
    % Convert unew to cell array format
    unewc = cell(1,options.n_dims);
    colons = repmat({':'},1,options.n_dims);
    for dimIndex = 1:options.n_dims
      unewc{dimIndex} = unew(colons{:},dimIndex);
    end
    % Scale the deformation to a unit-pixel grid
    uoldtmp = uoldc;
    if ~isempty(uoldtmp)
      for dimIndex = 1:options.n_dims
        % uoldtmp{dimIndex} = qinterp_grid_inverse(uoldtmp{dimIndex}/blocksz(dimIndex),'reflect');
        uoldtmp{dimIndex} = uoldtmp{dimIndex}/blocksz(dimIndex);
      end
    end
    unewtmp = unewc;
    for dimIndex = 1:options.n_dims
      unewtmp{dimIndex} = unewc{dimIndex}/blocksz(dimIndex);
    end
    % Evaluate the penalty
    if calc_grad
      [~,~,val,grad] = register_logdetpenalty_composition(uoldtmp,unewtmp,options.c);
    else
      [~,~,val] = register_logdetpenalty_composition(uoldtmp,unewtmp,options.c);
    end
    % Scale by lambda (and undo the unit-pixel scaling)
    val = options.lambda*val/prod(gridsz);
    if calc_grad
      for dimIndex = 1:options.n_dims
        grad{dimIndex} = (options.lambda/blocksz(dimIndex))*grad{dimIndex}/prod(gridsz);
      end
    end
  else
    % Skipping the regularization; we have to pre-allocate storage for the
    % output
    val = 0;
    if calc_grad
      grad = cell(1,options.n_dims);
      for dimIndex = 1:options.n_dims
        grad{dimIndex} = zeros(gridsz,options.class);
      end
    end
  end

  %% debugging code
  if ~isfinite(val)
    warning('val is NaN');
  end
%   val1 = val;
%   val1 / options.lambda
  
  %% The intensity-based distance/data penalty
  unew = reshape(unew,[n_grid options.n_dims]);
  for gridIndex = 1:n_grid
    center = ceil(size(mismatch{gridIndex})/2);
    if calc_grad
      % [thisval,thisgradc] = qinterp_grid(center+unew(gridIndex,:),mismatch{gridIndex},'nan');
      % thisgrad = cat(2,thisgradc{:});
      [thisval,thisgrad] = imqinterp(center+unew(gridIndex,:),mismatch{gridIndex});
    else
      % thisval = qinterp_grid(center+unew(gridIndex,:),mismatch{gridIndex},'nan');
      thisval = imqinterp(center+unew(gridIndex,:),mismatch{gridIndex});
    end
    % debugging code
    if ~isfinite(thisval)
      warning('thisval is NaN');
    end
    val = val+thisval;
    if calc_grad
      for dimIndex = 1:options.n_dims
        grad{dimIndex}(gridIndex) = grad{dimIndex}(gridIndex)+thisgrad(dimIndex);
      end
    end
  end
  
  %% debugging code
%   val2 = val - val1
  
  %% The edge barrier penalty
  if(barrier_mu > 0)
    BPmu = barrier_mu;
    for gridIndex = 1:n_grid    
      if calc_grad
        [BPsum, BPGrad] = register_barrierpenalty_composition(unew(gridIndex,:),mismatch{gridIndex});
      else
        [BPsum, ~] = register_barrierpenalty_composition(unew(gridIndex,:),mismatch{gridIndex});
      end
      %% debugging code
      if ~isfinite(BPsum)
        warning('BPsum is NaN');
      end
      val = val - BPmu * BPsum;
      if calc_grad
        for dimIndex = 1:options.n_dims
          grad{dimIndex}(gridIndex) = grad{dimIndex}(gridIndex) - ...
            BPmu * BPGrad(dimIndex);
        end
      end   
    end
    
    %% debugging code
%     val3 = val - val1 - val2;
%     val3 / BPmu
    
  end
  
  %% Convert to array format
  if calc_grad
    gradout = cat(options.n_dims+1,grad{:});
  end
  
  