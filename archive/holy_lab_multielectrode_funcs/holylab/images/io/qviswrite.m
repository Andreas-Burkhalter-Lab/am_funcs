function qviswrite(filename,data,options)
% QVISWRITE: save image in QVis format
%
% OpenQVis is a 3d volume rendering suite: http://openqvis.sourceforge.net/
% This file format is also the best-supported format for ImageVis3D.
%
% Syntax:
%   qviswrite(filename,data,options)
% where
%   filename is a string giving the name of the output header file (.dat
%     will be automatically appended if necessary)
%   data is a numeric array containing your image (can be of arbitrary
%     dimensionality), a stackmm object, or a string containing the name of
%     a raster data file (e.g., a .cam file). If you use the latter, you
%     must also supply the fields "precision" and "size" of options.
%   options is an optional structure which may have the following fields:
%     pixel_spacing: a vector of length equal to the # of spatial
%       dimensions in the image, specifying the physical displacement
%       corresponding to one grid point along each coordinate axis (for
%       example, a 3-d image sampled at 0.5um intervals in x,y and 5um
%       intervals along z would use [0.5 0.5 5] for this field)
%     precision: 'uint8' or 'uint16' are the two supported choices
%     scalefactor: values are multiplied by this value before converting to
%       a different precision
%     size: a 3-vector specifying the size
%     machineformat: a character passed to fwrite allowing you to change
%       the endian status (see fwrite)
%
% Note:  when you supply a stackmm object, if the setting for precision
% and/or machineformat is different from that of the native .cam file, the
% raw data are converted to the new format.
%
% Copyright 2010-2011 by Timothy E. Holy

  %% Check the inputs
  if (nargin < 3)
    options = struct;
  end
  [~,~,endian] = computer;
  options = default(options,...
    'pixel_spacing',[1 1 1],...
    'machineformat',lower(endian),...
    'scalefactor',1);
  if ischar(data)
    if ~isfield(options,'precision') || ~isfield(options,'size')
      error('If supplying a string for data, must supply the precision and size in options.')
    end
    sz = options.size;
  elseif isnumeric(data)
    options = default(options,'precision',class(data));
    sz = size(data);
  elseif isa(data,'stackmm')
    % stackmm object
    header = data.header;
    options = default(options,'precision',header.prec);
    rewrite = ~isequal(options.precision,header.prec || isempty(strmatch(lower(options.machineformat),{lower(endian),'n'},'exact')));
    sz = data.size;
    sz = sz(1:3); % fixme: write stack 1
  else
    error('data type not recognized');
  end

  if isempty(strmatch(options.precision,{'uint8','uint16'},'exact'))
    error('Precision must be uint8 or uint16');
  end
  
  sz(end+1:3) = 1;  % fill out to 3 spatial dimensions, required by qvis
  n_colors = 1;
  if (length(sz) > 3)
    n_colors = sz(4);
    if (n_colors ~=3 && n_colors ~=4)
      error('The number of colors must be 3 or 4');
    end
    sz = sz(1:3);
  end
  
  %% Parse the filename
  [pth,basename] = fileparts(filename);
  if ~isempty(pth)
    basename = [pth filesep basename];
  end
    
  %% Write the header file
  [fid,msg] = fopen([basename '.dat'],'w');
  if (fid < 0)
    error(msg);
  end
  if isnumeric(data)
    fprintf(fid,'ObjectFileName: %s.raw\n',basename);
  else
    fprintf(fid,'ObjectFileName: %s\n',data);
  end
  fprintf(fid,'TaggedFileName: -/home/tim/matlabfunc--\n');
  fprintf(fid,'Resolution: ');
  fprintf(fid,'%d ',sz(1:3));
  fprintf(fid,'\n');
  fprintf(fid,'SliceThickness: ');
  fprintf(fid,'%g ',options.pixel_spacing);
  fprintf(fid,'\n');
  if strcmp(options.precision,'uint8')
    if (n_colors > 1)
      fprintf(fid,'Format: UCHAR%d\n',n_colors);
    else
      fprintf(fid,'Format: UCHAR\n');
    end
  else
    if (n_colors > 1)
      fprintf(fid,'Format: USHORT%d\n',n_colors);
    else
      fprintf(fid,'Format: USHORT\n');
    end
  end
  fprintf(fid,'NbrTags: 0\n');
  fprintf(fid,'ObjectType: TEXTURE_VOLUME_OBJECT\n');
  fprintf(fid,'ObjectModel: RGBA\n');
  fprintf(fid,'GridType: EQUIDISTANT\n');
  fclose(fid);
  
  %% Write the data file
  if isnumeric(data)
    [fid,msg] = fopen([basename '.raw'],'w');
    if (fid < 0)
      error(msg);
    end
    if (n_colors > 1)
      data = permute(data,[4 1 2 3]);
    end
    fwrite(fid,data,class(data),0,options.machineformat);
    fclose(fid);
  elseif isa(data,'stackmm')
    if rewrite
      [fid,msg] = fopen([basename '.raw'],'w');
      if (fid < 0)
        error(msg);
      end
      for i = 1:sz(3)
        im = data(:,:,i,1);
        im = cast(options.scalefactor*im,options.precision);
        fwrite(fid,im,options.precision,0,options.machineformat);
      end
      fclose(fid);
    end
  end
end
