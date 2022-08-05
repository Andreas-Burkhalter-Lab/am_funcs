function gout = register_composeg(g,gcorr)
% REGISTER_COMPOSEG: Compose one deformation with another
% Given a deformation g and a second deformation gcorr, compute the
% combined deformation
%   gout(x) = g(gcorr(x)).
% 
% Syntax:
%   gout = register_composeg(g,gcorr)
%
% Implementation notes:
% This uses interpolation. Edge effects are handled by
% linear extrapolation.
%
% See also: REGISTER_G0.

% Copyright 2006-2007 by Timothy E. Holy
  
  inputs_numeric = isnumeric(g);
  if inputs_numeric
    g = g_array2cell(g);
  end
  if isnumeric(gcorr)
    gcorr = g_array2cell(gcorr);
  end
  n_dims = length(g);
  gout = cell(1,n_dims);
  szg = size(g{1});
  szg(end+1:n_dims) = 1;
  sz1flag = szg > 1;
  for dimIndex = 1:n_dims
    if (isempty(g{dimIndex}));
      gout{dimIndex} = gcorr{dimIndex};
    else
      if sz1flag(dimIndex)
        gout{dimIndex} = iminterp(g{dimIndex},gcorr{sz1flag},'extrap');
      else
        gout{dimIndex} = g{dimIndex} + gcorr{dimIndex} - 1;
      end
    end
  end
  if inputs_numeric
    gout = cat(n_dims+1,gout{:});
  end
