function savevolumevtk(filename, array, options)
%  savevolumevtk: Save an array in VTK format.
%  savevolumevtk(filename, array) saves a 3-D array of scalars in VTK format.
%  savevolumevtk(filename, array, options), where options may have the following fields:
%    pixel_spacing: a 3-vector specifying the spacing along each dimension
%      (default [1 1 1])
%    variable_name: a string giving the vtk name of the variable (default
%      'raw')
%
% The "scalars" can be alternatively be RGB-valued, in which case they need
%   to be uint8s in the shape [sz(1) sz(2) sz(3) 3].

% Copyright 2010 by Timothy E. Holy

  if (nargin < 3)
    options = struct;
  end
  options = default(options,'pixel_spacing',[1 1 1],'variable_name','raw');
  sz = size(array);
  rgb = (length(sz) > 3 && sz(4) == 3);
  [fid,msg] = fopen(filename, 'wt');
  if (fid < 0)
    error(msg);
  end
  fprintf(fid, '# vtk DataFile Version 2.0\n');
  fprintf(fid, 'Comment goes here\n');
  fprintf(fid, 'BINARY\n');
  fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
  fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', sz(1:3));
  fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
  fprintf(fid, 'SPACING    %f   %f   %f\n',options.pixel_spacing);
  fprintf(fid, '\n');
  fprintf(fid, 'POINT_DATA   %d\n', prod(sz(1:3)));
  if ~rgb
    fprintf(fid, 'SCALARS %s float\n',options.variable_name);
    fprintf(fid, 'LOOKUP_TABLE greys\n');
    fprintf(fid, '\n');
    fwrite(fid,array,'single');
  else
    fprintf(fid, 'COLOR_SCALARS %s 3\n',options.variable_name);
    fprintf(fid, '\n');
    fwrite(fid,permute(array,[4 1 2 3]),'uchar');
  end
  fclose(fid);
return