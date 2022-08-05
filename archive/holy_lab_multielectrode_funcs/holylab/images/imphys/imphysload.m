function ip = imphysload(fileinfo)
% IMPHYSLOAD: Virtually load image data from disk
% Syntax:
%   ip = imphysload(fileinfo)
% where
%   ip is the output IMPHYS structure---the image is
%     loaded into the image field of the structure;
% and
%   fileinfo contains the information about the file. This is a structure
%   which may have the following fields (see IMFILE for automated
%   ways of obtaining these):
%     fh  (a filehandle, see FILEHANDLE)
%     format (image format, e.g. 'raw', 'tif', 'jpeg', ...)
%     stacknum (the frame/stack number) 
%     stacktime (the time in secs of acquisition, relative to start)
%     size (read size of object)
%     prec (data type, e.g. 'uint16')
%     stimulus (valve #, if stimulus info is available)
%    + the following which are useful only for 'raw' images:
%     fpos (file position of image data in file)
%     pixelorder (Column-row ([1 2]) or row-column ([2 1])) ?
%
% See also: IMFILE, IMLOAD, FILEHANDLE.

% Copyright 2005 by Timothy E. Holy

  fileinfo = ipload_defaults(fileinfo);

  % Allocate all at once for speed
  ip = repmat(imphys,1,length(fileinfo));
  
  % Copy over relevant information
  for i = 1:length(fileinfo)
    ip(i).filename = fileinfo(i).fh.filename;
    ip(i).stacknum = fileinfo(i).stacknum;
  end
  if isfield(fileinfo,'stacktime')
    for i = 1:length(fileinfo)
      ip(i).stacktime = fileinfo(i).stacktime;
    end
  end
  if isfield(fileinfo,'stimulus')
    for i = 1:length(fileinfo)
      ip(i).stimulus = fileinfo(i).stimulus;
    end
  end

  % Set up virtual image structure
  vtmp = vimage(1,length(fileinfo));  % Pre-allocate
  for i = 1:length(fileinfo)
    ip(i).image = vtmp(i);
    push(ip(i).image,'imload',fileinfo(i));
  end
  
function fileinfo = ipload_defaults(fileinfo)
  if ~isfield(fileinfo,'abspathstr')
    [fileinfo.abspathstr] = deal('');
  end
  if ~isfield(fileinfo,'machfmt')
    [fileinfo.machfmt] = deal('n');
  end
  if ~isfield(fileinfo,'uselfs')
    [fileinfo.uselfs] = deal(1);
  end
  if ~isfield(fileinfo,'format')
    [fileinfo.format] = deal('');
  end
  if ~isfield(fileinfo,'fpos')
    [fileinfo.fpos] = deal([]);
  end
  if ~isfield(fileinfo,'prec')
    [fileinfo.prec] = deal('');
  end
  if ~isfield(fileinfo,'framenum')
    [fileinfo.framenum] = deal([]);
  end