function ucomp = register_phasecorr_composeu(uold,unew,chunksz)
% register_phasecorr_composeu: compose deformations for phasecorr registration
%
% Syntax:
%   ucomp = register_phasecorr_composeu(uold,unew,chunksz)
% where
%   uold is either a cell array or a numeric array specifying the old
%     deformation, in units of pixel shifts for the image chunks.
%   unew is a numeric array specifying the correction to the old
%     deformation
%   chunksz is the size of each image chunk, i.e., size(image)./szu (where
%     szu is the size of the spatial grid)
% and
%   ucomp is the composed deformation, expressed as a cell array.
%
% See also: register_phasecorr_improve.
  
% Copyright 2010 by Timothy E. Holy
  
  %% Parse inputs
  usz = size(unew);
  n_dims = usz(end);
  colons = repmat({':'},1,n_dims);
  if iscell(uold)
    uold = cat(n_dims+1,uold{:});
  end
  
  %% Handle unwarped inputs
  if isempty(uold)
    % We were starting from an unwarped image, so we just use the new one
    ucomp = cell(1,n_dims);
    for dimIndex = 1:n_dims
      ucomp{dimIndex} = unew(colons{:},dimIndex);
    end
    return
  end
  
  %% Compose the new and old
  chunksz = reshape(chunksz,[ones(1,n_dims) n_dims]);
  uold = bsxfun(@rdivide,uold,chunksz); % fractional shift
  unew = bsxfun(@rdivide,unew,chunksz);
  usz = size(unew);
  g0 = register_g0(usz(1:end-1),class(unew));
  gold = cell(1,n_dims);
  gnew = cell(1,n_dims);
  for dimIndex = 1:n_dims
    gold{dimIndex} = g0{dimIndex} + uold(colons{:},dimIndex);
    gnew{dimIndex} = g0{dimIndex} + unew(colons{:},dimIndex);
  end
  gcomp = register_composeg(gold,gnew);
  % Convert back to u representation (in pixel units)
  ucomp = cell(1,n_dims);
  for dimIndex = 1:n_dims
    ucomp{dimIndex} = (gcomp{dimIndex} - g0{dimIndex}) ...
      .* chunksz(dimIndex);
  end

  