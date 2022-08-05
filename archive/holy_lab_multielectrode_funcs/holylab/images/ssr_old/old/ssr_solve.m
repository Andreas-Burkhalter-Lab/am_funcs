function x = ssr_solve(H,b,features,options)
  stack_dims = length(features.stacksize);
  % The next two variables are to handle MATLAB's weirdness about 1d
  stack_dims_md = max(2,stack_dims);
  stack_dims_sd = sum(features.stacksize > 1);
  colons = {':'};
  colons = repmat(colons,1,stack_dims);
  n_features = length(features.S);
  n_pairs = size(features.overlappingPairsIndex,1);
  n_blocks = length(features.block);
  n_stacks = length(features.stacknum);
  %
  % Solve for the temporal components
  %
  if options.calcT
    % Note b.T is a nfeatures-by-nstacks matrix
    x.T = H.T \ b.T;
  end
  %
  % Solve for the spatial components
  %
  if options.calcS
    % Preallocation, with marking by NaNs to detect bad indexing
    x.S = features.S;
    for featureIndex = 1:n_features
      x.S{featureIndex} = nan(size(features.S{featureIndex}));
    end
    for blockIndex = 1:n_blocks
      % Split RHS data into blocks
      n_features_in_block = length(features.block(blockIndex).index);
      block_size = unique((features.block(blockIndex).features_coords_max - ...
        features.block(blockIndex).features_coords_min + 1)','rows');
      bblock = zeros([n_features_in_block block_size]);
      for featureBlockIndex = 1:n_features_in_block
        % Copy over the data in the given block
        for dimIndex = 1:stack_dims
          feature_coords{dimIndex} = features.block(blockIndex).features_coords_min(dimIndex,featureBlockIndex): ...
            features.block(blockIndex).features_coords_max(dimIndex,featureBlockIndex);
        end
        bblock(featureBlockIndex,colons{:}) = b.S{features.block(blockIndex).index(featureBlockIndex)}(feature_coords{:});
      end
      bblock = reshape(bblock,[n_features_in_block prod(block_size)]);
      % Conceptually we want to do this:
      %     xSblock = H.S{blockIndex} \ bblock;
      % but if the matrix is singular, we run into trouble.
      [tmpV,tmpD] = eig(H.S{blockIndex});
      Ddiag = diag(tmpD);
      isbad = (Ddiag < 1e-8*max(Ddiag));
      Ddiag(isbad) = Inf;
      tmpDinv = diag(1./Ddiag,0);
      xSblock = (tmpV*tmpDinv*tmpV') * bblock;
      xSblock = reshape(xSblock,[n_features_in_block block_size]);
      for featureBlockIndex = 1:n_features_in_block
        % Copy back the spatial profile into the block
        featureIndex = features.block(blockIndex).index(featureBlockIndex);
        for dimIndex = 1:stack_dims
          feature_coords{dimIndex} = features.block(blockIndex).features_coords_min(dimIndex,featureBlockIndex): ...
            features.block(blockIndex).features_coords_max(dimIndex,featureBlockIndex);
        end
        x.S{featureIndex}(feature_coords{:}) = xSblock(featureBlockIndex,colons{:}); %reshape?
      end
    end
  end
  %
  % Solve for the registration
  %
  % The calculations on particular stacks gets integrated into the
  % block-tridiagonal solver
  if options.calcRegistration
    if features.lambda
      % We are using regularization
      c = b.Registration;  % temporary storage "pre-allocation"
      xR = b.Registration; % output "pre-allocation"
      D = ssr_solve_calcdiag(H.registration,features,1);
      lamR = features.lambda * median(abs(H.registration.hessval_spatial)) * ...
        H.regularization;
      D = D - 2*lamR; % add in the diagonal component of the regularization
      Q = D;
      G = Q\lamR;
      G_all = cell(1,n_stacks);
      G_all{1} = G;
      offset0 = features.registration.offset(:,:,1);
      offsetp1 = features.registration.offset(:,:,1);
      rhs = b.Registration(:,1) - lamR*(offsetp1(:) - offset0(:));
      c(:,1) = Q\rhs;
      n_stacks = size(b.Registration,2);
      for stackIndex = 2:n_stacks
        D = ssr_solve_calcdiag(H.registration,features,stackIndex);
        D = D - 2*lamR; % Include the diagonal component of the regularization
        Q = D - lamR*G;
        G = Q\lamR;
        G_all{stackIndex} = G;
        offsetm1 = offset0;
        offset0 = offsetp1;
        if (stackIndex < n_stacks)
          offsetp1 = features.registration.offset(:,:,stackIndex+1);
          rhs = b.Registration(:,stackIndex) - lamR*(offsetp1(:)+offsetm1(:)-2*offset0(:));
        else
          rhs = b.Registration(:,stackIndex) - lamR*(offsetm1(:)-offset0(:));
        end
        c(:,stackIndex) = Q\(rhs - lamR*c(:,stackIndex-1));
      end
      xR(:,end) = c(:,end);
      for stackIndex = n_stacks-1:-1:1
        xR(:,stackIndex) = c(:,stackIndex) - G_all{stackIndex}*c(:,stackIndex+1);
      end
    else
      % We are not using regularization, so each stack is independent of
      % the others
      xR = nan(size(b.Registration));
      for stackIndex = 1:n_stacks
        D = ssr_solve_calcdiag(H.registration,features,stackIndex);
        xR(:,stackIndex) = D \ b.Registration(:,stackIndex);
      end
    end
    x.offset = reshape(xR,[stack_dims_sd n_features n_stacks]);
  end
  
  
function D = ssr_solve_calcdiag(registration,features,stackIndex)
  hessval = registration.hessval_spatial .* ...
	    features.T(registration.Tindex1,stackIndex) .* ...
	    features.T(registration.Tindex2,stackIndex);
  D = sparse(registration.rowIndex,registration.colIndex,...
	     hessval);