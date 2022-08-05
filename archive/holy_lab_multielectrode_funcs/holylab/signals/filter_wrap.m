function h = filter_wrap(h0,sz,center)
% filter_wrap: put a multidimensional filter into position for zero-shift
% Syntax:
%   hout = filter_wrap(hin,sz)
%   hout = filter_wrap(hin,sz,center)
% where
%   hin is the input filter
%   sz is a vector containing the sizes of the dimensions of the object(s)
%     to be filtered
%   center is the position of the peak/center of the filter (default: in
%     the center of the supplied region)
% and
%   hout is the wrapped filter, of size sz, in which the center of the
%   filter has been moved to the corner and the other terms wrapped using
%   periodic boundary conditions.

% Copyright 2008 by Timothy E. Holy

  nd = ndims(h0);
  szh = size(h0);
  szh0 = szh;
  if (length(sz) < length(szh))
    % Assume one-dimensional
    szh = szh(szh > 1);
    nd = length(szh);
  end
  if (nargin < 3)
    center = floor(szh/2);
  end
  if (nd > 1)
    h = zeros(sz,class(h0));
  else
    szh1 = szh0;
    szh1(szh0 > 1) = sz;
    h = zeros(szh1,class(h0));
  end
  x = cell(1,nd);
  % Create a wrap-around filter, so that it's centered on the center
  for i = 1:nd
    x{i} = [(1:center(i))+sz(i)-center(i) (center(i)+1:szh(i))-center(i)];
  end
  h(x{:}) = h0;
  return


  
  %% Old code (useless?)
%   szu = sz(sz>1); % sz with unity-dimensions suppressed
%   ndims = length(szu);
%   psfsz = size(psfi);
%   psfszu = psfsz(psfsz > 1);
%   if (length(psfszu) ~= ndims)
%     error('The source and target have different dimensionality');
%   end
% 
%   % Re-shape the PSF so that the center is mapped to the corner; this
%   % will keep filtering operations from introducing a global shift.
%   psf = zeros(sz,class(psfi));
%   half = cell(ndims,2);
%   halfsize = zeros(1,ndims);
%   for i = 1:ndims
%     halfsize(i) = floor(psfszu(i)/2);
%     half{i,1} = 1:halfsize(i);
%     half{i,2} = halfsize(i)+1:psfszu(i);
%   end
%   from_coord = cell(1,ndims);
%   to_coord = cell(1,ndims);
%   for i = 0:2^ndims-1
%     for j = 1:ndims
%       b = bitget(i,j); % 0 = first half, 1 = 2nd half
%       from_coord{j} = half{j,b+1};
%       to_coord{j} = (1-b)*szu(j) - halfsize(j) + from_coord{j};
%     end
%     psf(to_coord{:}) = psfi(from_coord{:});
%   end
