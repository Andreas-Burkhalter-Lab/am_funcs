function nrrdwrite(filename,data,options)
% NRRDWRITE: save image in NRRD format
%
% NRRD, or "nearly-raw raster data," is a format used by the Utah
% Scientific Imaging and Visualization group. This format is imported by
% tools such as ImageVis3d.
%
% Syntax:
%   nrrdwrite(filename,data,options)
% where
%   filename is a string giving the name of the file (.nrrd will be
%     automatically appended if necessary)
%   data is a numeric array containing your image (can be of arbitrary
%     dimensionality)
%   options is an optional structure which may have the following fields:
%     pixel_spacing: a vector of length equal to the # of spatial
%       dimensions in the image, specifying the physical displacement
%       corresponding to one grid point along each coordinate axis (for
%       example, a 3-d image sampled at 0.5um intervals in x,y and 5um
%       intervals along z would use [0.5 0.5 5] for this field)
%     kind: '' or 'scalar' indicates that data is to be interpreted as a
%       scalar intensity, 'rgb' or 'rgba' specifies that the last dimension
%       corresponds to color information.

% Copyright 2010 by Timothy E. Holy

  %% Parse inputs
  if (nargin < 3)
    options = struct;
  end
  options = default(options,'kind','');
  
  %% Parse the filename and open the file for writing
  [pth,basename] = fileparts(filename);
  if isempty(pth)
    outfile = [basename '.nrrd'];
  else
    outfile = [pth filesep basename '.nrrd'];
  end
  
  [fid,msg] = fopen(outfile,'w');
  if (fid < 0)
    error(msg);
  end
  closeFile = onCleanup(@() fclose(fid));  % so the file gets closed at exit, no matter what
  
  %% Write the header
  % Specify the data format
  fprintf(fid,'NRRD0004\n');
  switch class(data)
    case 'single'
      fprintf(fid,'type: float\nendian: little\n');
    case 'double'
      fprintf(fid,'type: double\nendian: little\n');
    case 'uint8'
      fprintf(fid,'type: uint8\n');
    case 'uint16'
      fprintf(fid,'type: uint16\n');
  end
  
  % Specify the data size
  sz = size(data);
  sz = sz(sz > 1);
  n_dims = length(sz);
  fprintf(fid,'dimension: %d\n',n_dims);
  fprintf(fid,'sizes: ');
  fprintf(fid,'%d ',sz);
  fprintf(fid,'\n');
  
  % Write tags that interpret the dimensions (the "kind" tag)
  switch(lower(options.kind))
    case {'','scalar'}
      n_spatial_dims = n_dims;
      kind = '';
    case 'rgb'
      n_spatial_dims = n_dims-1;
      kind = [options.kind '-color'];
      if (sz(end) ~= 3)
        error('If specifying RGB, the last dimension must have size 3');
      end
    case 'rgba'
      n_spatial_dims = n_dims-1;
      kind = [options.kind '-color'];
      if (sz(end) ~= 3)
        error('If specifying RGB, the last dimension must have size 3');
      end
  end
  fprintf(fid,'kind: ');
  for i = 1:n_spatial_dims
    fprintf(fid,'space ');
  end
  fprintf(fid,'%s\n',kind);

  % Specify the pixel spacing
  if isfield(options,'pixel_spacing')
    fprintf(fid,'space directions: ');
    if (length(options.pixel_spacing) ~= n_spatial_dims)
      error('format:dimensions','Mismatch in the number of spatial dimensions and the pixel_spacing.\nDid you forget to specify RGB or RGBA?');
    end
    for i = 1:n_spatial_dims
      fprintf(fid,'(');
      for j = 1:i-1
        fprintf(fid,'0,');
      end
      fprintf(fid,'%g',options.pixel_spacing(i));
      for j = i+1:n_spatial_dims
        fprintf(fid,',0');
      end
      fprintf(fid,') ');
    end
    fprintf(fid,'\n');
  end
  
  fprintf(fid,'encoding: raw\n\n');
  
  %% Write the raw data
  fwrite(fid,data,class(data));

  % closing is done by the cleanup function above
end
