function rng = nanbox(A)
% nanbox: find the largest subarray containing no NaNs
%
% This algorithm is designed for image registration, where the NaNs arise
% from interpolation and are therefore localized to the edges. It is
% equivalent to successively cutting off slices (the one containing the
% most NaNs) from the edges of the array until there are no more NaNs.
%
% Syntax:
%   rng = nanbox(A)
% where
%   A is an array which may contain some NaNs (each NaN should be connected
%     to the edge of the array via a "path" of NaNs)
% and
%   rng is a n_dims-by-2 matrix indicating the coordinate ranges of the
%     largest box with no NaN.
%
% Example:
%   A = [nan(1,4) 4; nan 1 2 3 4; nan nan 2 3 4; nan 1 2 3 nan];
% Then
%   rng = nanbox(A);
% returns  [2 3; 3 5].
% You can get the subregion this way:
%   for i = 1:size(rng,1)
%     rngc{i} = rng(i,1):rng(i,2);
%   end
%   Asnip = A(rngc{:});

% Copyright 2010 by Timothy E. Holy

  n_dims = ndims(A);
  sz = size(A);
  colons = repmat({':'},1,n_dims);
  rng = zeros(n_dims,2);
  s = cell(1,n_dims);  % will hold count of NaNs along each coordinate
  for dimIndex = 1:n_dims
    rng(dimIndex,:) = [1 sz(dimIndex)];
    s{dimIndex} = countnans(A,dimIndex);
  end
  
  isdone = false;
  while ~isdone
    % We'll try to keep the box as big as possible, so compute its size
    rngsz = diff(rng,1,2)+1;
    % Find the edge that has the most NaNs
    dim = 1;
    sedge = s{1}(rng(1,:));
    [smax,index] = max(sedge);
    sztop = rngsz(1);
    for dimIndex = 2:n_dims
      sedge = s{dimIndex}(rng(dimIndex,:));
      [smaxtmp,indextmp] = max(sedge);
      if (smaxtmp > smax || (smaxtmp == smax && rngsz(dimIndex) > sztop))
        dim = dimIndex;
        index = indextmp;
        smax = smaxtmp;
        sztop = rngsz(dimIndex);
      end
    end
    % Count the # of NaNs eliminated by clipping this slice, and subtract
    % from the total count
    sliceIndex = rng(dim,index);
    colonstmp = colons;
    colonstmp{dim} = sliceIndex;
    for dimIndex = 1:n_dims
      sslice = countnans(A(colonstmp{:}),dimIndex);
      if (dimIndex == dim)
        s{dimIndex}(sliceIndex) = 0;
      else
        s{dimIndex} = s{dimIndex} - sslice;
      end
    end
    A(colonstmp{:}) = 0;  % clear the NaNs we just counted
    % Clip this slice from the range
    if (index == 1)
      rng(dim,index) = rng(dim,index)+1;
    else
      rng(dim,index) = rng(dim,index)-1;
    end
    % Test to see if we are done
    isdone = true;
    for dimIndex = 1:n_dims
      isdone = ~any(s{dimIndex}>0);
      if ~isdone
        break
      end
    end
  end
  