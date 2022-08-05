function [eigenvec,eigenval] = lda(X,group,options)
% LDA: Linear Discriminant Analysis
%
% Syntax:
%   [eigenvec,eigenval] = lda(X,group)
% where
%   X is a d-by-N matrix, where N is the number of data points (each of
%     which is d-dimensional); thus, each data point is a column of X
%   group is a 1-by-N vector, containing a number that indicates the
%     identity of each data point.  Most typically, these will be
%     integers ranging from 1 to n_groups.
%   eigenvec is d-by-(n_groups-1)
%   eigenval is 1-by-(n_groups-1)
%
% Alternative syntax:
%   [eigenvec,eigenval] = lda(X,means)
% where
%   X is as above
%   means is d-by-m (m = # of cluster centers), means should be shifted by
%     the data average. Note that means does not have to be derived from
%     the same samples used in X.
%
% Alternatively,
%   [eigenvec,eigenval] = lda(...,options)
% permits more control. The option fields are:
%   maxpergroup (default Inf): if set to something finite, a subset of the
%     supplied data will be used to perform LDA, choosing at most
%     maxpergroup points in each group. (This syntax requires the
%     lda(X,group) syntax.)
%   regularized (default true): if true, the eigenvalue calculation will
%     be done to avoid useless null directions (this is expecially
%     useful/important for small sample sizes, i.e., N < d or when N is of
%     order d). Our approach follows Yu & Yang's "Direct LDA" (2001).
%
% See also: PCA.

% Copyright 2006-2007 by Timothy E. Holy

  if (nargin < 3)
    options = struct;
  end
  options = default(options,'regularized',true,'maxpergroup',Inf);
  
  [d,N] = size(X);
  if (numel(group) == N)
    % X, group syntax
    % Convert the group labels to 1,2,3,...
    [ugroup,tmp,group_label] = unique(group);
    [clabel,nlabel] = agglabel(group_label);
    % Check to see if we are supposed to be using only a subset of the
    % points
    for groupIndex = 1:length(nlabel)
      if (nlabel(groupIndex) > options.maxpergroup)
        skip = ceil(nlabel(groupIndex)/options.maxpergroup);
        clabel{groupIndex} = clabel{groupIndex}(1:skip:end);
        nlabel(groupIndex) = length(clabel{groupIndex});
      end
    end
    % Calculate the group means
    n_groups = length(nlabel);
    mu = zeros(d,n_groups);
    for groupIndex = 1:n_groups
      mu(:,groupIndex) = sum(X(:,clabel{groupIndex}),2);
    end
    xbar = sum(mu,2)/sum(nlabel);
    mu = mu ./ repmat(nlabel,d,1);
    % Calculate the within-group displacement from the mean
    dX = cell(1,n_groups);
    for groupIndex = 1:n_groups
      dX{groupIndex} = X(:,clabel{groupIndex}) - ...
        repmat(mu(:,groupIndex),1,nlabel(groupIndex));
    end
    dX = cat(2,dX{:});
    % Finally, subtract off the overall mean
    dmu = mu - repmat(xbar,1,n_groups);
  else
    % X, means syntax
    xbar = mean(X,2);
    if (size(group,1) ~= d)
      error('Can''t determine which syntax is being used, check sizes of inputs');
    end
    dX = X - repmat(xbar,1,N);
    dmu = group;
    options.regularized = true;
  end
  
  if options.regularized
    % Roundoff is a real problem if we straightforwardly try to maximize
    % SB/Stot. So we instead diagonalize SB and then use those dimensions
    % to diagonalize SW.
    % Begin by diagonalizing the between-groups covariance
    [UB,SB] = svd(dmu,'econ');
    sB = diag(SB);
    ikeep = sB >= 1e-8*max(sB);  % throw out near-zero singular values
    Y = UB(:,ikeep);
    Z = Y*diag(1./sB(ikeep));
    % Project the within-groups (or total) "half-covariance" (covariance is
    % dX*dX') using Z, to restrict the dimensions to those that span the
    % means
    projW = Z'*dX;
    [UW,SW] = svd(projW,'econ');
    sW = diag(SW);
    eigenval = 1./sW.^2;
    eigenvec = Z*UW*diag(1./sW);
    % Start with largest eigenvalues first
    eigenval = eigenval(end:-1:1);
    eigenvec = eigenvec(:,end:-1:1);
    eigenvec = eigenvec ./ repmat(sqrt(sum(eigenvec.^2,1)),d,1);
  else
    SB = zeros(d,d);
    SW = zeros(d,d);
    for groupIndex = 1:n_groups
      dmean = mu(:,groupIndex) - xbar;
      SB = SB + nlabel(groupIndex) * (dmean*dmean');
      Xtmp = X(:,clabel{groupIndex}) - repmat(mu(:,groupIndex),[1 ...
        nlabel(groupIndex)]);
      SW = SW + Xtmp * Xtmp';
    end
    [eigenvec,eigenval] = eig(SB,SW);
    eigenval = diag(eigenval);
    [eigenval,sort_order] = sort(eigenval,'descend');
    eigenvec = eigenvec(:,sort_order);
    eigenval(n_groups:end)=[];
    eigenvec(:, n_groups:end)=[];
  end
  