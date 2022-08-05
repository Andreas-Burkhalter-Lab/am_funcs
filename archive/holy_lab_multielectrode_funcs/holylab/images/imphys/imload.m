function im = imload(fileinfo)
% IMLOAD: Load real image data from disk
% This is the "real" cousin to IMPHYSLOAD, which virtually loads the
% images (sets up a vimage).
%
% Syntax:
%   im = imload(fileinfo)
% where
%   im is the image
% and
%   fileinfo contains the information about the file(s). This is a structure
%   which should have the following fields (see IMFILE for automated
%   ways of obtaining these):
%     fh (a filehandle class object, see FILEHANDLE).
%     format (image format, e.g. 'raw', 'tif', 'jpeg', ...)
%     single (if true, a single image/file)
%     stacknum (the frame/stack number) 
%     stacktime (the time in secs of acquisition, relative to start)
%     size (read size of object)
%     prec (data type, e.g. 'uint16')
%     stimulus (valve #, if stimulus info is available)
%    + the following which are useful only for 'raw' images:
%     fpos (file position of image data in file)
%     pixelorder (Column-row ([1 2]) or row-column ([2 1])) ?
%
% See also: IMPHYSLOAD, IMFILE, FILEHANDLE.

% Copyright 2005 by Timothy E. Holy

  if (length(fileinfo) > 1)
    error('May only load one image at a time');
  end

  filename = [fileinfo.fh.abspathstr fileinfo.fh.filename];
  % Load raw format
  if strcmp(fileinfo.format,'raw')
    imnpix = prod(fileinfo.size);
    [fh,msg] = fopen(fileinfo.fh,'r');
    if (fh.fid < 0)
      error(sprintf(['Can''t open file %s. Check permissions? ' ...
                     'Message:\n%s'],filename,msg));
    end
    if fh.uselfs
      if strcmp(fileinfo.prec,'uint16')
        im = readuint16lfs(fh.fid,imnpix,[0 0],fileinfo.fpos);
      else
        error(['Raw precision ' fileinfo.prec ' not implemented']);
      end
    else
      status = fseek(fh.fid,fileinfo.fpos,'bof');
      if (status < 0)
          error('Bad seek, dude');
      end
      [im,count] = fread(fh.fid,imnpix,['*' fileinfo.prec]);
      if (count < prod(imnpix))
          error('Too few values were read from the file; file truncated?');
      end
    end
    fclose(fh)
    % Reshape data to yield something of the right size and conventional
    % orientation
    if ~isempty(fileinfo.pixelorder)
      im = reshape(im,fileinfo.size(fileinfo.pixelorder));
      if ~issorted(fileinfo.pixelorder)
        im = permute(im,fileinfo.pixelorder);
      end
    end
    
    % special handling for the raw data created by Andor software
    %ttDataFileName=fileinfo.fh.filename;
    %ttHeader=imreadheader(imdatafilename2headerfilename(ttDataFileName));
    if(strcmp(fileinfo.camera, 'Andor 434') && strcmp(fileinfo.format,'raw'))
       im=im(end:-1:1,:); % reverse row order
    end
    
    return
  end
  % Read Metamorph stk format
  if strcmp(fileinfo.format,'stk')
    data = tiffread(filename, fileinfo.stackpos);
    im = data.image;
    return
  end
  % Read other formats
  if fileinfo.single
    % For single image filetypes
    im = imread(filename,fileinfo.format);
  else
    % For multi-tiffs, etc.
    im = imread(filename,fileinfo.format,fileinfo.stacknum);
  end
  