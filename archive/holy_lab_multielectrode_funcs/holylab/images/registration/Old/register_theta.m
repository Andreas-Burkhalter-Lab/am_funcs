function theta = register_theta(imM)
% REGISTER_THETA: calculate the grid mask for the moving image
% Syntax:
%   theta = register_theta(imM)

  imMtmp = squeeze(imM);
  n_dims = ndims(imMtmp);
  % Find NaNs and their nearest-neighbors
  theta = ~isnan(imMtmp);
  theta = imerode(theta,ones(repmat(3,1,n_dims)));
  % Fill in the edges
  theta = killedges(theta,0);
  % Do it again to give a buffer zone of width 2
  theta = imerode(theta,ones(repmat(3,1,n_dims)));
  theta = cast(theta,class(imM));
  theta = reshape(theta,size(imM));
  