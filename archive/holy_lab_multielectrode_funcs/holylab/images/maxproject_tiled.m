function im = maxproject_tiled(tmm,tilerng,dim)
% maxproject_tiled: perform max projection with tiled representations
%
% Syntax:
%   im = maxproject_tiled(tmm,tilerng,dim)
% where
%   tmm is a tiledmm object;
%   tilerng is a cell array, {coord1min:coord1max,coord2min:coord2max,...},
%     with 4 elements
%   dim is the dimension along which to calculate the max projection.
%
% On output, im is an array containing the max projection along dimension
% dim for the specified set of tiles.
%
% See also: cell2mat.

% Copyright 2012 by Timothy E. Holy

  %% Parse args
  if ~iscell(tilerng) || length(tilerng) ~= 4
    error('tilerng must be a cell array with 4 elements');
  end
  if ~isscalar(dim)
    error('dim must be an integer');
  end
  
  ntiles = tmm.ntiles;
  for i = 1:4
    if isequal(tilerng{i},':')
      tilerng{i} = 1:ntiles(i);
    end
  end
  tileinddim = tilerng{dim};
  tilerng{dim} = 1;
  tileind = cell(1,4);
  [tileind{:}] = ndgrid(tilerng{:});
  ntiles = size(tileind{1});
  ntiles(end+1:4) = 1;
  N = numel(tileind{1});
  
  %% Calculate maxproject across tiles
  imt = cell(ntiles);
  first = true;
  tic;
  t = 0;
  for i = 1:N
    thisind = {tileind{1}(i),tileind{2}(i),tileind{3}(i),tileind{4}(i)};
    fullind = thisind;
    fullind{dim} = tileinddim;
    tmp = tmm.tiles(fullind{:});
    tmp = max(tmp,[],ndims(tmp));
    imt{thisind{:}} = max(tmp,[],dim);
    if (toc-t > 5)
      if first
        fprintf('Progress (%% done): ');
        first = false;
      end
      fprintf('%d...',round(100*i/N));
      t = toc;
    end
  end
  if ~first
    fprintf('done.\n');
  end
  
  %% Combine results from tiles
  im = cell2mat(imt);
  % Truncate any extra pixels due to tiling
  sz = tmm.size;
  imsz = size(im);
  imsz(end+1:4) = 1;
  mx = min(imsz,sz);
  for i = 1:4
    thisind{i} = 1:mx(i);
  end
  im = im(thisind{:});
end
