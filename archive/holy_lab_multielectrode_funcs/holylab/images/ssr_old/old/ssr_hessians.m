function H = ssr_hessians(features,options)
  stack_dims = length(features.stacksize);
  % The next two variables are to handle MATLAB's weirdness about 1d
  stack_dims_md = max(2,stack_dims);
  stack_dims_sd = sum(features.stacksize > 1);
  n_features = length(features.S);
  n_pairs = size(features.overlappingPairsIndex,1);
  n_blocks = length(features.block);

  if options.calcT
    %
    % The Hessian for the temporal optimization
    %
    hessval_od = nan(n_pairs,1);
    for pairIndex = 1:n_pairs
      tmpIndex = features.overlappingPairsIndex(pairIndex,:);
      tmpRange = features.overlappingPairsCoords(pairIndex,:);
      for dimIndex = 1:stack_dims
        coords1{dimIndex} = tmpRange{1}(dimIndex,1):tmpRange{1}(dimIndex,2);
        coords2{dimIndex} = tmpRange{2}(dimIndex,1):tmpRange{2}(dimIndex,2);
      end
      tmpProductImage = features.S{tmpIndex(1)}(coords1{:}) .* ...
        features.S{tmpIndex(2)}(coords2{:});
      hessval_od(pairIndex) = sum(tmpProductImage(:));
    end
    hessval_d = nan(n_features,1);
    for featureIndex = 1:n_features
      tmpProductImage = features.S{featureIndex}.^2;
      hessval_d(featureIndex) = sum(tmpProductImage(:));
    end
    rowIndex = [features.overlappingPairsIndex(:,1); ...
      features.overlappingPairsIndex(:,2); ...
		  (1:n_features)'];
    colIndex = [features.overlappingPairsIndex(:,2); ...
		  features.overlappingPairsIndex(:,1); ...
		  (1:n_features)'];
    hessval = [hessval_od; hessval_od; hessval_d];
    H.T = sparse(rowIndex,colIndex,hessval,n_features,n_features);
  end
  if options.calcS
    %
    % The Hessian for the spatial optimization
    %
    % Here there isn't just one Hessian, there is one for each block of
    % pixels with common feature membership.
    % This uses the zero-neighbor-drift approximation (meaning that
    % overlapping features "don't" drift apart from each other)
    for blockIndex = 1:n_blocks
      n_features_in_block = length(features.block(blockIndex).index);
      H.S{blockIndex} = nan(n_features_in_block,n_features_in_block);
      for featureBlockIndex1 = 1:n_features_in_block
        featureIndex1 = features.block(blockIndex).index(featureBlockIndex1);
        for featureBlockIndex2 = featureBlockIndex1+1:n_features_in_block
          featureIndex2 = features.block(blockIndex).index(featureBlockIndex2);
          H.S{blockIndex}(featureBlockIndex1,featureBlockIndex2) = ...
            sum(features.T(featureIndex1,:) .* features.T(featureIndex2,:));
          H.S{blockIndex}(featureBlockIndex2,featureBlockIndex1) = ...
            H.S{blockIndex}(featureBlockIndex1,featureBlockIndex2);
        end
        H.S{blockIndex}(featureBlockIndex1,featureBlockIndex1) = ...
          sum(features.T(featureIndex1,:).^2);
      end
    end
  end
  if options.calcRegistration
    %
    % The Hessian for the registration
    %
    % This does not include any of the regularization terms; nor does it
    % factor in the temporal terms, which will be integrated during the
    % inversion process (for reasons of memory conservation).
    dS = ssr_feature_gradients(features.S,stack_dims_sd);
    hessval_od = nan(stack_dims_sd,stack_dims_sd,n_pairs);
    for pairIndex = 1:n_pairs
      tmpIndex = features.overlappingPairsIndex(pairIndex,:);
      tmpRange = features.overlappingPairsCoords(pairIndex,:);
      for dimIndex = 1:stack_dims
        % Compute the overlapping ranges
        coords1{dimIndex} = tmpRange{1}(dimIndex,1):tmpRange{1}(dimIndex,2);
        coords2{dimIndex} = tmpRange{2}(dimIndex,1):tmpRange{2}(dimIndex,2);
      end
      for dimIndex1 = 1:stack_dims_sd
        for dimIndex2 = 1:stack_dims_sd
          tmpProductImage = dS{dimIndex1,tmpIndex(1)}(coords1{:}) .* ...
            dS{dimIndex2,tmpIndex(2)}(coords2{:});
          hessval_od(dimIndex1,dimIndex2,pairIndex) = sum(tmpProductImage(:));
        end
      end
    end
    hessval_d = nan(stack_dims_sd,stack_dims_sd,n_features);
    for featureIndex = 1:n_features
      for dimIndex1 = 1:stack_dims_sd
        for dimIndex2 = dimIndex1:stack_dims_sd
          tmpProductImage = dS{dimIndex1,featureIndex} .* ...
            dS{dimIndex2,featureIndex};
          hessval_d(dimIndex1,dimIndex2,featureIndex) = sum(tmpProductImage(:));
          hessval_d(dimIndex2,dimIndex1,featureIndex) = ...
            hessval_d(dimIndex1,dimIndex2,featureIndex);
        end
      end
    end
    
    % Now comes the indexing to create the (sparse) hessian. We have to
    % convert the (dim,dim,pairIndex) (for off-diagonal) and
    % (dim,dim,featureIndex) (for diagonal) into a single index.
    % The off-diagonal (pairs) contribution:
    indexP1=repmat(features.overlappingPairsIndex(:,1)',[stack_dims_sd^2 1]);
    indexP2=repmat(features.overlappingPairsIndex(:,2)',[stack_dims_sd^2 1]);
    indexDimPairs1 = repmat((1:stack_dims_sd),[stack_dims_sd n_pairs]);
    indexDimPairs2 = repmat((1:stack_dims_sd)',[1 stack_dims_sd*n_pairs]);
    indexP1 = indexP1(:);
    indexP2 = indexP2(:);
    indexDimPairs1 = indexDimPairs1(:);
    indexDimPairs2 = indexDimPairs2(:);
    singleIndexP1 = sub2ind([stack_dims_sd n_features],indexDimPairs1,indexP1);
    singleIndexP2 = sub2ind([stack_dims_sd n_features],indexDimPairs2,indexP2);
    % The "diagonal" contribution (self-pairs)
    indexD =repmat(1:n_features,[stack_dims_sd^2 1]);
    indexD = indexD(:);
    indexDimFeatures1 = repmat((1:stack_dims_sd),[stack_dims_sd n_features]);
    indexDimFeatures2 = repmat((1:stack_dims_sd)',[1 stack_dims_sd*n_features]);
    indexDimFeatures1 = indexDimFeatures1(:);
    indexDimFeatures2 = indexDimFeatures2(:);
    singleIndexD1 = sub2ind([stack_dims_sd n_features], indexDimFeatures1,indexD);
    singleIndexD2 = sub2ind([stack_dims_sd n_features], indexDimFeatures2,indexD);
    hessval_od_transpose = permute(hessval_od,[2 1 3]);
    
    rowIndex = [singleIndexP1; singleIndexP2; singleIndexD1];
    colIndex = [singleIndexP2; singleIndexP1; singleIndexD2];
    hessval = [hessval_od(:); hessval_od_transpose(:); hessval_d(:)];
    
    % Don't go all the way to the sparse matrix, because the temporal
    % factors will be integrated later. But do set up the indexing to make
    % that easy.
    H.registration.rowIndex = rowIndex;
    H.registration.colIndex = colIndex;
    H.registration.hessval_spatial = hessval;
    H.registration.Tindex1 = [indexP1; indexP2; indexD];
    H.registration.Tindex2 = [indexP2; indexP1; indexD];
  end    
  if options.calcRegularization
    % This computes the "interesting parts" of the Hessian for the
    % following penalty:
    %    E = (1/2) lambda * sum_t sum_(i,j nbrs) [(b_i(t+1)-b_i(t)) -
    %                                             (b_j(t+1)-b_j(t))]^2
    % This is the only temporally-off-diagonal block in the Hessian, as
    % it operates between adjacent stacks.  However, it also contains
    % diagonal elements
    % Note each i,j block is diagonal in the dimensions; hence we don't
    % repmat over stack_dims^2 like we do in the registration Hessian
    indexP1=repmat(features.overlappingPairsIndex(:,1)',[stack_dims_sd 1]);
    indexP2=repmat(features.overlappingPairsIndex(:,2)',[stack_dims_sd 1]);
    indexDim = repmat((1:stack_dims_sd)',[1 n_pairs]);
    indexP1 = indexP1(:);
    indexP2 = indexP2(:);
    indexDim = indexDim(:);
    singleIndexP1 = sub2ind([stack_dims_sd n_features],indexDim,indexP1);
    singleIndexP2 = sub2ind([stack_dims_sd n_features],indexDim,indexP2);
    H.regularization = sparse([singleIndexP1; singleIndexP2],...
      [singleIndexP2; singleIndexP1],...
      2*ones(2*length(singleIndexP1),1),...
      stack_dims_sd*n_features,stack_dims_sd*n_features,...
      2*length(singleIndexP1)+stack_dims_sd*n_features);
    H.regularization = spdiags(-sum(H.regularization,2),0,H.regularization);
    % This term represents the elementary component of the
    % regularization. Each row of blocks looks like this:
    %  [0 0 0 ... 0 R -2*R R 0 ... 0 0 0] 
    % where R = lambda * H.regularization
  end