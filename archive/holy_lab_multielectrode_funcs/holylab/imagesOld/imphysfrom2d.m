function ip = imphysfrom2d(imfile,headerfile)
% IMPHYSFROM2D: create an imphys structure from a file with 2d images
% Syntax:
%   ip = imphysfrom2d(imfile,headerfile)
%
% See also: IMPHYS.

% Copyright 2004 by Timothy E. Holy

   [pathstr,basename,extname] = fileparts(imfile);
   if (nargin < 2 || isempty(headerfile)) && strcmp(extname,'.tif')
      % todo: some info, such as spf, is missing from tif file
      
      info=imfinfo(imfile);
      for i=1:length(info)
         ip(i).imfile=imfile;
         ip(i).imfilefmt = 'tif';
         ip(i).imfilemachfmt = '';
         ip(i).imfileprec = '';
         ip(i).camera='unknown';
         ip(i).headerfile = '';
         ip(i).pixelorder = '';
         ip(i).date = info(i).FileModDate;
         ip(i).height = info(i).Height;
         ip(i).width = info(i).Width;
         ip(i).depth = 1;
         ip(i).xrange = [1 info(i).Width];
         ip(i).yrange = [1 info(i).Height];
         ip(i).zrange = [1 1];
         ip(i).stacknum = i;
         ip(i).z2um = NaN;
         ip(i).stimulus = [];
         ip(i).stacktime = [];
         ip(i).stackfpos = [];
         ip(i).imrange=[info(i).MinSampleValue info(i).MaxSampleValue]; % todo: or from .BitDepth?
      end
      return;
   end
   
  h = imreadheader(headerfile);
  
  % todo: set imrange field
  
  if(strcmp(extname,'.tif'))
     h.filefmt = 'tif'; 
     info=imfinfo(imfile);
     h.nstacks=length(info);
  else
     % treate it as raw
     tFileSize=filesize(imfile);
     h.nstacks=floor(tFileSize/h.nbytes/h.height/h.width);
  end
  
  % if(isempty(extname))
  
  % Convert stimulus information into a lookup
  stim = [];
  if ~isempty(h.stim)
    stim = h.stim(1,1);
    for i = 1:(size(h.stim,1)-1)
      stim = [stim zeros(1,diff(h.stim([i i+1],2)))+h.stim(i,1)];
    end
    stim = [stim zeros(1,h.nstacks-length(stim))+h.stim(end,1)];
  end
  % Fill in fields
  for i = 1:h.nstacks
    ip(i).imfile = imfile;
    ip(i).imfilefmt = h.filefmt;
    ip(i).imfilemachfmt = h.machfmt;
    ip(i).imfileprec = h.prec;
    ip(i).camera=h.camera;
    ip(i).headerfile = headerfile;
    ip(i).pixelorder = h.pixelorder;
    ip(i).date = h.date;
    ip(i).height = h.height;
    ip(i).width = h.width;
    ip(i).depth = 1;
    ip(i).xrange = [1 h.width];
    ip(i).yrange = [1 h.height];
    ip(i).zrange = [1 1];
    ip(i).stacknum = i;
    ip(i).z2um = NaN;
    if ~isempty(stim)
      ip(i).stimulus = stim(i);
    end
    ip(i).stacktime = h.stacktime(i);
    ip(i).stackfpos = h.nbytes*h.width*h.height*(i-1);
  end
