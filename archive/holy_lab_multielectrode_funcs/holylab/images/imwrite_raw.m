function [fileinfo,nextfpos] = imwrite(fid,image,mode)
% IMWRITE: save an image in raw format
% Syntax:
%   fileinfo = imwrite(fid,image)
% where
%   fid is the file identifier
%   image is image, either real or virtual
% If image is virtual (a vimage), then after this call returns you might
% want to consider adding a new method for loading the data from this
% file. You do this by executing a push(image,'imload',fileinfo)
% after this function returns.
%
% Alternative syntax:
%   [fileinfo,nextfpos] = imwrite(fpos,image,'fake')
% does a 'fake' write, advancing the "file position" (just a counter) to
% its next position, dependent upon the number of bytes needed to write
% image.
%
% See also: IMLOAD.

% Copyright 2005 by Timothy E. Holy

  if (nargin < 3)
    mode = '';
  end
  
  fileinfo = struct('format','raw','pixelorder',[1 2]);
  if strcmp(mode,'fake')
    fileinfo.fpos = fid;
    if isa(image,'vimage')
      fileinfo.prec = imclass(image);
    else
      fileinfo.prec = class(image);
    end
    fileinfo.size = size(image);
    nextfpos = fid + sizeof(fileinfo.prec)*prod(fileinfo.size);
    % Note filename and machfmt need to be filled in by calling function
    return
  end
  [fileinfo.fh,fileinfo.fpos] = filehandle(fid);
  if isa(image,'vimage')
    tmp = eval(image);
  else
    tmp = image;
  end
  fileinfo.size = size(tmp);
  prec = class(tmp);
  count = fwrite(fid,tmp,prec);
  if (count < numel(tmp))
    error('Write failed!');
  end
  fileinfo.prec = prec;
