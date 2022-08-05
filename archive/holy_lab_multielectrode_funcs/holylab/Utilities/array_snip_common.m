function varargout = array_snip_common(varargin)
  % array_snip_common: compute coordinates for snipping out regions of arrays
  %
  % In tasks like image registration, one frequently needs to pull out a
  % block from the image, and do so in a way that the excised block has a
  % consistent relationship with a block (perhaps of a different size)
  % pulled out from a different image. This function calculates
  % common "destination" coordinates from the input coordinates.
  %
  % Syntax:
  %   [spanout,sz] = array_snip_common(spanin)
  % where
  %   spanin is a cell array, one for each dimension. Each element is a
  %     n-by-2 matrix, where we plan to make n cutouts. Each row of this
  %     matrix contains [coordmin coordmax] for this dimension.
  % and
  %   spanout is a cell array with similar organization, giving the
  %     coordinates 
  %
  % Example:
  %   im = imread('cameraman.tif');
  %   c1src = {20:120,80:180};  % Coords used for the first snippet
  %   c2src = {65:130,125:170}; % Coords used for the second snippet
  %   [c1dest,c2dest,sz] = array_snip_common(c1src,c2src);
  %   im1 = nan(sz);            % Create NaN-padded first snippet
  %   im2 = nan(sz);            % Create NaN-padded second snippet
  %   im1(c1dest{:}) = im(c1src{:});   % Stuff the snippet in place
  %   im2(c2dest{:}) = im(c2src{:});
  %   figure; imshowrgb(im1,im2);
  %
  % See also: coords2spans, spans2coords.
  
  % Copyright 2012 by Timothy E. Holy
  
  n_arrays = length(varargin);
  if nargout < n_arrays
    error('You need at least as many outputs as inputs (perhaps one extra for sz)');
  end
  
  [minsrc,maxsrc] = coords2spans(varargin{:});  % convert to an array representation
  offset = min(minsrc,[],1)-1;
  mindest = bsxfun(@minus,minsrc,offset);
  maxdest = bsxfun(@minus,maxsrc,offset);
  cdest = cell(1,n_arrays);
  [cdest{:}] = spans2coords(mindest,maxdest);
  varargout(1:n_arrays) = cdest;
  if nargout > n_arrays
    varargout{n_arrays+1} = max(maxdest,[],1);
  end
end
