function fileinfo = imfile(imfile,headerfile)
% IMFILE: prepare to read image data from a file
% This function doesn't actually load the image data, it simply sets up
% the file information needed to read it, for both "standard" and lab
% image formats.
%
% Syntax:
%   ip = imfile(imfile)
% where imfile is the filename for a multi-tif
%   ip = imfile(imfile,headerfile)
% where imfile is the filename multi-raw or multi-tif
% and headerfile is the name of the text file describing the images
%   ip = imfile(imfiles)
% where imfiles is a cell array of single image filenames, with any
% standard image format.
%   ip = imfile(imfiles,fileinfo)
% is like the above, except you can specify 'raw' image files (using this
% extension) and provide information about the dimensions, etc. in
% the structure fileinfo.  The fields of this structure are described in
% IMPHYSLOAD; note that often only the dimensions and precision will need
% to be specified. The input fileinfo may either be a structure (in which
% case the choices apply to all) or a structure array (if you need to
% customize).
%
% See also: IMPHYSLOAD, IMLOAD.

% Copyright 2005 by Timothy E. Holy

  % Case: cell array of single images
  if iscell(imfile)
    nfiles = length(imfile);
    fileinfo(nfiles).stacknum = 0;  % Pre-allocate for speed
    for i = 1:nfiles
      if ischar(imfile{i})
        [pathstr,basename,extname] = ipf2_abspath(imfile{i});
        fh = filehandle('abspathstr',pathstr,...
          'filename',[basename extname]);
        fileinfo(i).fh = fh;
        fileinfo(i).stacknum = i;
        if strcmp(extname,'.raw')
          fileinfo(i).format = 'raw';
          fileinfo(i).fpos = 0;
          % User had better supply size, prec, etc!
          if (nargin > 1)
            j = i;
            if (length(headerfile) == 1)
              j = 1;
            end
            fileinfo(i).size = headerfile(j).size;
            fileinfo(i).prec = headerfile(j).prec;
            fileinfo(i).pixelorder = headerfile(j).pixelorder;
          end
        else
          info = imfinfo(imfile{i});
          fileinfo(i).format = info.Format;
          fileinfo(i).size = [info.Height info.Width];
        end
        fileinfo(i).single = 1;
      else
        error('Syntax not recognized');
      end
    end
    return;
  end

  [pathstr,basename,extname] = ipf2_abspath(imfile);
  
  % Case: multi-tif, no custom header info
  % Metamorph stk format might take a bit of doing? Or does imfinfo work here?
  if (nargin < 2 || isempty(headerfile)) && ...
        strcmp(extname,'.tif')
    % Note: some info, such as spf, is missing from tif file
    info=imfinfo(imfile);
    fileinfo(length(imfile)).stacknum = 0;  % Pre-allocate for speed
    fh = filehandle('abspathstr',pathstr,...
                    'filename',[basename extname]);
    for i=1:length(info)
      fileinfo(i).fh = fh;
      fileinfo(i).stacknum = i;
      fileinfo(i).format = 'tif';
      fileinfo(i).size = [info(i).Height info(i).Width];
      fileinfo(i).single = 0;
    end
    return;
  end
  
  % OK, we have one of our custom header files available
  h = imreadheaderOld(headerfile);
  
  % todo: set imrange field
  
  if(strcmp(extname,'.tif'))
    h.filefmt = 'tif'; 
    info=imfinfo(imfile);
    h.nstacks=length(info);
  elseif ~isfield(h,'nstacks')
    % treate it as raw
    tFileSize=filesize(imfile);
    h.nstacks=floor(tFileSize/h.nbytes/h.height/h.width);
  end
  
  % Convert stimulus information into a lookup
  stim = [];
  if ~isempty(h.stim)
    %stim = h.stim(1,1);
    for i = 1:(size(h.stim,1)-1)
      stim = [stim zeros(1,diff(h.stim([i i+1],2)))+h.stim(i,1)];
    end
    stim = [stim zeros(1,h.nstacks-length(stim))+h.stim(end,1)];
  end

  % Fill in fields
  fileinfo(h.nstacks).stacknum = 0;      % Pre-allocate for speed
  uselfs = 1;
  if ispc
    uselfs = 0;
  end
  fh = filehandle('abspathstr',pathstr,...
                  'filename',[basename extname],...
                  'machfmt',h.machfmt,...
                  'uselfs',uselfs);
  for i = 1:h.nstacks
    fileinfo(i).fh = fh;
    fileinfo(i).format = h.filefmt;
    fileinfo(i).stacknum = i;
    fileinfo(i).stacktime = h.stacktime(i);
    fileinfo(i).prec = h.prec;
    if isfield(h,'depth')
      fileinfo(i).size = [h.width h.height h.depth];
    else
      fileinfo(i).size = [h.width h.height];
    end
    fileinfo(i).pixelorder = h.pixelorder;
    fileinfo(i).camera = h.camera;
    if strcmp(h.filefmt,'raw')
      fileinfo(i).fpos = h.nbytes*prod(fileinfo(i).size)*(i-1);
    end
    fileinfo(i).single = 0;
    if ~isempty(stim)
      fileinfo(i).stimulus = stim(i);
    end
  end

  
function [pathstr,basename,extname] = ipf2_abspath(filename)
% Get absolute path
  if isunix
    if strcmp(filename(1),filesep)
      % It's already the absolute path
      [pathstr,basename,extname] = fileparts(filename);
    else
      % Specified as relative path, get absolute
      [pathstr,basename,extname] = fileparts([pwd filesep filename]);
    end
    pathstr = [pathstr filesep];
  else
    warning('Windows not fully implemented')
    % What to do on Windows?
    [pathstr,basename,extname] = fileparts(filename);
  end
  
