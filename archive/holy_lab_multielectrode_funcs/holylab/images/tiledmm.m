classdef tiledmm
%  tiledmm:  memory-mapped access to tiled arrays
%
% This class is similar to stackmm, but with on-disk organization of the
% data. While stackmm accesses data that are stored in "array" format
% (column-major ordering), tiledmm accesses data stored as a
% "tiled array," meaning that the large array has been split up into a
% set of non-overlapping cubes. For accessing data in a small region of
% interest, this representation can have significant performance
% advantages.
%
%  Syntax:
%    tmm = tiledmm(basename)
%    tmm = tiledmm(basename, 'bias', b)
%  where
%    basename is a string containing the name of the tiled array (with
%      extension ".tiled") and its header file (extenstion ".tiledhdr").
%
% You access data as if tmm is a large array, just as with stackmm:
%    roi = tmm(25:55,200:325,3:8,:);
% A tiledmm object will likely be more efficient than a stackmm object if
% you are trying to get information about a small local volume,
% particularly if you are interested in what happens there over time. It
% will likely be less efficient than a stackmm if you're trying to get a
% whole stack, or frame, at a time.
%
% You can also extract the raw tile data in this fashion:
%    imt = tmm.tiles(1:5,1:6,1:5,1:16);
% imt will be an 8-dimensional array, consisting of the tiles indexed in
% parentheses. Or, if you want the whole tiles corresponding to a
% coordinate range of the data, you can say
%    ims = tmm.astiles(25:55,200:325,3:8,:);
% This will return a strcture containing the information needed to convert
% this to an array (the tile data, plus the array offsets for the
% "snipped-out" region). 
% These two methods of accessing the data are faster than asking for an
% array.
%
% In addition to the same operations supported by stackmm, you also have:
%    tmm.tilesize: the size of each tile
%    tmm.ntiles: the number of tiles along each coordinate dimension
% You cannot change the setting of these fields, as they are determined by the
% layout on-disk.
%
% Unlike stackmm, you cannot chain multiple files together.
%
% See also: make_tiled_movie, stackmm.

%  Copyright 2012 by Timothy E. Holy
    
  properties (GetAccess = public, SetAccess = private)
    sz = [];
    tilesize = [];
    header = struct;
    filename = '';
    type = '';
    bias = [];
    badframes = [];
    badpixels = [];
  end
  properties (Access = private)
    mm = [];     % the memory-mapped object(s)
    n_stacks = [];
  end
  methods
    %% Constructor
    function tmm = tiledmm(basename, varargin)
      if isa(basename, 'tiledmm')
        % Copy Constructor
        tmm = basename;
      else
        % Create a new tiledmm
        if length(basename) > 9 && strcmp(basename(end-8:end),'.tiledhdr')
          basename = basename(1:end-9);
        end
        tmm.filename = basename;
        % Read the header and parse the fields
        hdrfilename = [basename '.tiledhdr'];
        if ismat(hdrfilename)
          s = load('-mat',hdrfilename);
          tmm.header = s.header;
          tmm.sz = s.arraysize;
          tmm.tilesize = s.tilesize;
          tmm.type = tmm.header.prec;
        else
          [fid,msg] = fopen(hdrfilename);
          if (fid < 0)
            error(msg)
          end
          s = fread(fid,'*char');
          tmm.sz = str2num(key2value(s,'arraysize')); %#ok<ST2NM>
          tmm.tilesize = str2num(key2value(s,'tilesize')); %#ok<ST2NM>
          try
            hdrfile = key2value(s,'header file');
            tmm.header = imreadheader(hdrfile);
            tmm.type = tmm.header.prec;
          catch
            tmm.type = key2value(s,'type');
          end
        end
        % Check that the file is the right size
        tiledfilename = [basename '.tiled'];
        dirtmp = dir(tiledfilename);
        n_tiles = ceil(tmm.sz ./ tmm.tilesize);
        n_bytes_per_tile = prod(tmm.tilesize)*sizeof(tmm.header.prec);
        if prod(n_tiles)*n_bytes_per_tile ~=  dirtmp.bytes
          error([tiledfilename ' is not of the expected size'])
        end
        % Set up the memory-mapped file
        tmm.mm = memmapfile(tiledfilename,...
              'format',{tmm.type,[tmm.tilesize n_tiles],'im'},...
              'repeat',1);
      end
    end
    function im = subsref(tmm,s)
      switch s(1).type
        case '()'
          %% Extract data
          % Parse the indices
          [c,l] = parse_indices(tmm,s(1).subs);
          % Convert to a linear index on the tile array
          n_tiles = ceil(tmm.sz ./ tmm.tilesize);
          tcoord = acoord2tilecoord(c,tmm.tilesize,n_tiles);
          % Grab the data
          im = reshape(tmm.mm.data.im(tcoord),l);
%           % This old version might be more efficient in a faster language!
%           % Pre-calculate the coordinate intersections
%           cmin = cellfun(@min,c);
%           cmax = cellfun(@max,c);
%           tilemin = ceil(cmin./tmm.tilesize);
%           tilemax = ceil(cmax./tmm.tilesize);
%           tile_coord = cell(1,n_dims);
%             
%           % Iterate over tiles
%           tilecur = tilemin;
%           tsnip = cell(1,n_dims);
%           asnip = cell(1,n_dims);
%           while tilecur(end) <= tilemax(end)
%             % For this tile, find the intersection of requested indices and
%             % tile indices
%             for idim = 1:n_dims
%               thistile = (tilecur(idim)-1)*tmm.tilesize(idim)+1:tilecur(idim)*tmm.tilesize(idim);
%               [thisarray,ic] = intersect(c{idim},thistile);
%               tsnip{idim} = thisarray - thistile(1) + 1;
%               asnip{idim} = ic;
%             end
%             if any(isempty(tsnip))
%               continue
%             end
%             tsnip_cat = cat(2,tsnip,num2cell(tilecur));
%             tmp = tmm.mm.data.im(tsnip_cat{:});
%             % Correct for camera bias
%             if ~isempty(tmp) && ~isempty(tmm.bias)
%               if isscalar(tmm.bias)
%                 tmp = tmp-tmm.bias;
%               else
%                 tmp = bsxfun(@minus,tmm.bias(tsnip{1:ndims(tmm.bias)}));
%               end
%             end
%             im(asnip{:}) = tmp;
%             % Increment the tile iterator
%             tilecur(1) = tilecur(1)+1;
%             idim = 1;
%             while idim < n_dims && tilecur(idim) > tilemax(idim)
%               tilecur(idim) = tilemin(idim);
%               idim = idim+1;
%               tilecur(idim) = tilecur(idim)+1;
%             end
%           end
        case '.'
          %% Other "get" methods that use subsref
          % These are necessary due to a limitation in Matlab's class
          % implementation; if you supply a subsref method so you can
          % overload parentheses (which we do for getting image data), then
          % you are also required to handle all other possible subsref
          % calls on your own.
          switch s(1).subs
            case 'tiles'
              if ~strcmp(s(2).type,'()')
                error('Syntax should be tmm.tiles(i,j,k,l)')
              end
              n_dims = length(tmm.sz);
              sub = cat(2,{':',':',':',':'},s(2).subs);
              im = tmm.mm.data.im(sub{:});
            case 'astiles'
              c = parse_indices(tmm,s(2).subs);
              n_dims = length(c);
              cmin = cellfun(@min,c);
              cmax = cellfun(@max,c);
              tilemin = ceil(cmin./tmm.tilesize);
              tilemax = ceil(cmax./tmm.tilesize);
              offset = (tilemin-1).*tmm.tilesize;
              coords = cell(1,n_dims);
              for i = 1:n_dims
                coords{i} = tilemin(i):tilemax(i);
              end
              coords = cat(2,repmat({':'},1,n_dims),coords);
              im = struct('start',cmin - offset,...
                'stop',cmax - offset,...
                'tiles',tmm.mm.data.im(coords{:}));
            case 'header'
              im = tmm.header;
            case 'size'
              im = tmm.sz;
            case 'tilesize'
              im = tmm.tilesize;
            case 'ntiles'
              im = ceil(tmm.sz./tmm.tilesize);
            case 'type'
              im = tmm.type;
            case 'filename'
              im = tmm.filename;
            case 'bias'
              im = tmm.bias;
            case 'badframes'
              im = tmm.badframes;
            case 'badpixels'
              im = tmm.badpixels;
            otherwise
              error(['Field ' s.subs ' not recognized']);
          end
      end
    end
    %% Assignment methods
    % These support the user setting "badframes" and "badpixels" after
    % construction
    function obj = subsasgn(obj,s,b)
      switch s.type
        case '.'
          switch s.subs
            case 'badframes'
              if ~isequal(size(b),obj.sz(3:4))
                error('The dimensions of the badframes matrix does not match');
              end
              obj.badframes = logical(b);
            case 'badpixels'
              if ~isequal(size(b),obj.sz(1:2))
                error('The dimensions of the badpixels matrix does not match');
              end
              obj.badpixels = logical(b);
            otherwise
              error(['Field ' s.subs ' not recognized']);
          end
        otherwise
          error(['subsasgn: ' s.type ' not recognized']);
      end
    end
  end
end

function [c,l] = parse_indices(tmm,subs)
  n_dims = length(tmm.sz);
  c = cell(1,n_dims);
  l = zeros(1,n_dims);
  for i = 1:n_dims
    if strcmp(subs{i},':')
      c{i} = 1:tmm.sz(i);
      l(i) = tmm.sz(i);
    else
      c{i} = subs{i};
      l(i) = length(c{i});
      if any(c{i} > tmm.sz(i))
        error('Tiledmm:outofrange',...
          'Array dimensions are only [%s]',num2str(tmm.sz));
      end
    end
  end
end

% Convert array coordinates to a linear index on the tiled representation
function tindex = acoord2tilecoord(ca,tilesize,nTiles)
  n_dims = length(ca);
  s = cumprod([1 tilesize nTiles]);
  s = s(1:end-1);
  tindex = mod(ca{1}(:)-1,tilesize(1))*s(1) + (ceil(ca{1}(:)/tilesize(1))-1)*s(1+n_dims);
  o = ones(1,n_dims);
  for i = 2:n_dims
    sz = o;
    sz(i) = length(ca{i});
    tmp = reshape(ca{i},sz);
    tindex = bsxfun(@plus,tindex,mod(tmp-1,tilesize(i))*s(i) + (ceil(tmp/tilesize(i))-1)*s(i+n_dims));
  end
  tindex = tindex+1;
end
