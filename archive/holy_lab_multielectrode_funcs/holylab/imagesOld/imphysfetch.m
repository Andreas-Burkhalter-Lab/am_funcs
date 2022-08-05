function im = imphysfetch(ip)
% IMPHYSFETCH: Return image data
% Syntax:
%   im = imphysfetch(ip)
% where
%   ip is the IMPHYS structure (must be a single element);
% and
%   im is the image.
%
% If the image is already present in the IMPHYS structure, the image is
% simply returned. Otherwise, the image is loaded from disk. Image
% transforms and cropping are performed before returning the image.
%
% See also: IMPHYS, IMPHYSFROM2D.
  
% Notes:
%   Consider making this more efficient by doing cropping by
%   altered file position and frame size data, never reading the unused
%   columns. This might necessitate giving up the tform support, or
%   at least checking to see if tform is present.
  
% Copyright 2004 by Timothy E. Holy
  need2crop = 1;
  uselfs = 0;
  if (length(ip) > 1)
    error('Will return only a single frame/stack')
  end
  if (isfield(ip,'image') & ~isempty(ip.image))
    im = ip.image;       % Return the pre-computed image
    return;
  end
  % OK, must load from disk
  if strcmp(ip.imfilefmt,'raw')
    % Determine size of the image
    imsize = [ip.height ip.width ip.depth];
    imnpix = prod(imsize);
    % Support a variety of data formats
    machfmt = 'n';
    if isfield(ip,'imfilemachfmt')
      machfmt = ip.imfilemachfmt;
    end
    prec = 'uint16';
    if isfield(ip,'imfileprec')
      prec = ip.imfileprec;
    end
    if (strcmp(machfmt,'n') & strcmp(prec,'uint16') & uselfs)
      % Use lfs versions to support files > 2GB
      fid = openlfs(ip.imfile);
      if (fid < 0)
        error(['Can''t open file ' ip.imfile]);
      end
      im = readuint16lfs(fid,imnpix,[0 0],ip.stackfpos);
      closelfs(fid);
    else
      % Use matlab versions for general support
      [fid,message] = fopen(ip.imfile,'r',machfmt);
      if (fid < 0)
        error(['Can''t open file ' ip.imfile ', message: ' message]);
      end
      fseek(fid,ip.stackfpos,'bof');
      im = fread(fid,imnpix,['*' prec]);
      fclose(fid);
    end
    % Reshape data to yield something of the right size and conventional
    % orientation
    im = reshape(im,imsize(ip.pixelorder));
    if ~issorted(ip.pixelorder)
      im = permute(im,ip.pixelorder);
    end
    % special handling for the raw data created by Andor software
    if(strcmp(ip.camera, 'Andor 434'))
       im=im(end:-1:1,:); % reverse row order
    end
  elseif strcmp(ip.imfilefmt,'tif')
    im = imread(ip.imfile,ip.imfilefmt,ip.stacknum,'PixelRegion',{ip.yrange,ip.xrange});
    need2crop = 0;
  end
  % If needed, transform the image
  if (isfield(ip,'tform') & ~isempty(ip.tform))
    im = imtransform(im,ip.tform,'XData',[1 ip.width],'YData',...
                     [1 ip.height]);
  end
  % If needed, crop the image
  if need2crop
     fn = {'xrange','yrange','zrange'};
     fn2 = {'width','height','depth'};
     for i = 1:3
        r{i} = ':';
        crng = ip.(fn{i});
        if any(crng ~= [1 ip.(fn2{i})])
           r{i} = crng(1):crng(2);
        end
     end
     if ( ~iscellstr(r) || length(strmatch(':',r,'exact')) < length(r) )
        if (ndims(im) == 2)
           im = im(r{[2 1]});
        else
           im = im(r{[2 1 3]}); % todo: this need to be verified
        end
     end % if need crop (r has double OR r has non-: )
  end