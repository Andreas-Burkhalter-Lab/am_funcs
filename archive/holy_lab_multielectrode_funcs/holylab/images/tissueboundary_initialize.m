function zdec = tissueboundary_initialize(stk,thresh,gridsep)
% TISSUEBOUNDARY_INITIALIZE: initial guess for tissue boundary
% Syntax:
%   z = tissueboundary_initialize(stk,thresh)
%   zdec = tissueboundary_initialize(stk,thresh,gridsep)
% where
%   stk is an image stack, with z as the _first_ coordinate, and up
%     corresponds to increasing z;
%   thresh is the threshold for the difference between tissue and black
%     space 
%   gridsep (optional) is the number of pixels between grid points; if
%     omitted, it returns the height for each (x,y) value
% and
%   z/zdec is the height as a function of x and y. zdec is the "decimated"
%     version, smoothed over the values for each x & y, that you get when
%     you supply something for gridsep.
%
% Note that the first syntax allows you to directly visualize the guessed
% surface using TISSUEBOUNDARY_CHECK.
%
% See also: TISSUEBOUNDARY, TISSUEBOUNDARY_CHECK.

% Copyright 2006 by Timothy E. Holy

  sz_stk = size(stk);
  sz_z = sz_stk(2:end);
  if isscalar(sz_z)
    sz_z = [1 sz_z];
  end
  z = ones(sz_z);
  for i = 1:numel(z)
      tmp = find(stk(:,i) > thresh,1,'last');
      if ~isempty(tmp)
          z(i) = tmp;
      end
  end
  if (nargin < 3)
    zdec = z;
    return
  end
  % Now decimate
  ngrid = ceil(sz_z./gridsep);
  for i = 1:length(sz_z)
    x{i} = round(linspace(1,sz_z(i),ngrid(i)+1));
    for j = 1:ngrid(i)
      X{i}{j} = x{i}(j):x{i}(j+1);
    end
  end
  % The next code is specific for d=1 or d=2
  if isvector(z)
    for i = 1:length(X{1})
      zblock = z(X{1}{i});
      zdec(i) = max(zblock);
    end
  else
    for i = 1:length(X{1})
      for j = 1:length(X{2})
        zblock = z(X{1}{i},X{2}{j});
        zdec(i,j) = max(zblock(:));
      end
    end
  end
