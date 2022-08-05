function make_tiled_movie(smm,options)
% make_tiled_movie: rearrange data in tiled format
%
% A tiled data format makes certain data-access operations much more
% efficient (and others, less efficient). See tiledmm for details.
%
% Syntax:
%   make_tiled_movie(smm)
%   make_tiled_movie(smm,options)
% where
%   smm is the stackmm object for the original file
%   options is a structure which may have the following fields:
%     tilesize (default [32 32 4 8]): the size of individual tiles in the
%       output file. Your choice for this setting has important
%       memory-usage implications: with the default, you must have enough
%       RAM to load 8 stacks at once. (The last non-singleton setting
%       determines the memory requirement, so [32 32 4 1] would only
%       require that 4 _frames_ be loaded at once.)
%       For reasons of performance during later reads, it's best to make
%       prod(tilesize)*sizeof(datatype) be a multiple of the page size,
%       4096 bytes. Beyond this requirement, there is no need to make it a
%       power of 2.
%     basename is the base name (i.e., without extension) of the output
%       file. If you do not supply this input, the file will have the same
%       base name as the smm input, with a '.tiled' extension for the data
%       and a '.tiledhdr' for the header.
%
% Tiled data can be conveniently accessed using tiledmm.
%
% See also: tiledmm.
  
% Copyright 2010-2011 by Timothy E. Holy
  
  %% Parse the inputs
  header = smm.header;
  if (nargin < 2)
    options = struct;
  end
  [~,basename] = fileparts(smm.filename);
  options = default(options,'basename',basename,'tilesize',[32 32 4 8]);
  options.tilesize(end+1:4) = 1;  % pad with 1s if necessary
  arraysize = smm.size;
  tilesize = options.tilesize;
  
  %% Write the header
  save([options.basename '.tiledhdr'],'header','arraysize','tilesize');
  
  %% Prepare the output data file
  [fid,msg] = fopen([options.basename '.tiled'],'w');
  if (fid < 0)
    error(msg);
  end
  cleanup_fid = onCleanup(@() fclose(fid));
  
  %% Determine whether we need to read multiple stacks per chunk
  if tilesize(4) == 1
    nChunks = [ceil(arraysize(3)/tilesize(3)) arraysize(4)];
    dimlookupChunk = [3 4];
    iSnipChunk1 = {':',':'};
    iSnipChunk2 = cell(1,2);
    subtilesize = tilesize(1:3);
    nTiles = ceil(arraysize(1:3)./subtilesize);
  else
    nChunks = ceil(arraysize(4)/tilesize(4));
    dimlookupChunk = 4;
    iSnipChunk1 = {':',':',':'};
    iSnipChunk2 = {[]};
    subtilesize = tilesize;
    nTiles = ceil(arraysize./tilesize);
  end
  nTiles(end) = 1; % per chunk, not for the whole experiment
  
  %% Iterate over chunks
  tic;
  tlast = toc;
  first = true;
  for iChunk = 1:prod(nChunks)
    % Prepare the indices used to snip out the chunk
    subChunk = ind2sub_matrix(nChunks,iChunk);
    for iDimChunk = 1:length(nChunks)
      iDim = dimlookupChunk(iDimChunk);
      iSnipChunk2{iDimChunk} = (subChunk(iDimChunk)-1)*tilesize(iDim)+1:min(subChunk(iDimChunk)*tilesize(iDim),arraysize(iDim));
    end
    iSnipChunk = cat(2,iSnipChunk1,iSnipChunk2);
    % Read the chunk
    stk = smm(iSnipChunk{:});
    % Write out the tiles
    split_chunk(nTiles,subtilesize);
    % Generate output for the user
    if (toc - tlast > 5)
      if first
        fprintf('%% done: ');
        first = false;
      end
      fprintf('%d%%...',round(100*iChunk/prod(nChunks)));
      tlast = toc;
    end
  end
  if ~first
    fprintf('done.\n');
  end
 
  function split_chunk(nTiles,tilesize)
    nDims = length(nTiles);
    cTilesize = num2cell(tilesize);  % for padding (when needed)
    nPerTile = prod(tilesize);
    indexSnip = cell(1,nDims);
    for iTile = 1:prod(nTiles)
      % Calculate the indices for this tile
      subTile = ind2sub_matrix(nTiles,iTile);
      for iDim = 1:nDims
        indexSnip{iDim} = (subTile(iDim)-1)*tilesize(iDim)+1:min(subTile(iDim)*tilesize(iDim),size(stk,iDim));
      end
      tile = stk(indexSnip{:});
      % If a tile extends beyond the original boundaries, pad with zeros
      if numel(tile) < nPerTile
        tile(cTilesize{:}) = 0;
      end
      % Save to disk
      count = fwrite(fid,tile,header.prec);
      if count < nPerTile
        error('Error writing tile. Could the disk be out of space?');
      end
    end
  end    
end
  