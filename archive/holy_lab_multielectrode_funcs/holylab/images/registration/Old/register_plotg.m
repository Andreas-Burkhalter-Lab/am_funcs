function register_plotg(g,n_vertices)
% register_plotg: display a deformation
% Syntax:
%   register_plotg(g)
%   register_plotg(g,n_vertices)
% where
%   g is a cell array defining a deformation, of the type you get back
%     from REGISTER_MULTIRESOLUTION or REGISTER_G0 (see the latter for a
%     detailed explanation of the format of g).
%   n_vertices (default 10) specifies the maximum number of vertices
%     along each coordinate.
%
% The plot has a marker at each grid point x, and a line connecting to
%   the "warped location" g(x).
%
% Obviously, the display for multidimensional deformations is
%   problematic; you can choose the dimensions that interest you by
%   supplying empty matrices for the dimensions you plan to leave out of
%   the plot (all but 2).  It selects the first element in all
%   "neglected" dimensions.
% You can, of course, manually create a g that includes whatever portion
%   you want to examine.
% TODO: allow this to create 3d plots
% 
% See also: REGISTER_G0, REGISTER_MULTIRESOLUTION.
  
% Copyright 2006 by Timothy E. Holy
  
  if (nargin < 2)
    n_vertices = 10;
  end
  n_dims = length(g);
  if (n_dims > 2)
    for dimIndex = 1:n_dims
      is_empty(dimIndex) = isempty(g{dimIndex});
    end
    coord_index = find(~is_empty,2,'first');
    colons = repmat({':'},1,n_dims);
    g = g(coord_index); % First two non-empty coordinates
    subset_flag = true(size(is_empty));
    subset_flag(coord_index) = false;
    colons(subset_flag) = {1};
    for dimIndex = 1:2
      g{dimIndex} = g{dimIndex}(colons{:});
    end
  end
  x = {1:size(g{1},1),1:size(g{1},2)};
  [X,Y] = ndgrid(x{:});
  if (ischar(n_vertices) && strcmp(n_vertices,'all'))
    indx1 = x{1};
    indx2 = x{2};
  else
    indx1 = unique(round(linspace(1,x{1}(end),n_vertices)));
    indx2 = unique(round(linspace(1,x{2}(end),n_vertices)));
  end
  %dg{1} = g{1} - X;
  %dg{2} = g{2} - Y;
  X = X(indx1,indx2);
  Y = Y(indx1,indx2);
  gX = g{1}(indx1,indx2);
  gY = g{2}(indx1,indx2);
  X = X(:);
  Y = Y(:);
  gX = gX(:);
  gY = gY(:);
  figure
  plot(Y,X,'ko')
  line([Y'; gY'],[X'; gX'],'Color','k')
  %quiverc(Y(indx1,indx2),X(indx1,indx2),dg{2}(indx1,indx2),dg{1}(indx1,indx2));
  set(gca,'YDir','reverse')
  axis image