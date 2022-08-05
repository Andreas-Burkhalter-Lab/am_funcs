function create_iv3d_script(scriptname,uvf_filename,pngbasename,rotatesequence,translatesequence,options)
% create_iv3d_script: write an ImageVis3d rendering script
%
% Syntax:
%   create_iv3d_script(scriptname,uvf_filename,pngbasename,rotatesequence,translatesequence,options)
% where
%   scriptname is the filename of the IV3D script (the output)
%   uvf_filename is the name of the .uvf file you want to load and render.
%     Set to empty if you want to load the data file manually (using the
%     debug window "execute" command to launch the script, e.g.
%       execute scriptname
%   pngbasename is the base name to use in writing the image frames, e.g.,
%     'frame'. The actual frames will be written to files such as
%     frame017.png.
%   rotatesequence is an n_frames-by-3 matrix, where each row specifies the
%     rotations to apply before rendering the next frame
%   translatesequence is an n_frames-by-3 matrix, each row specifying the
%     translation to be applied before rendering. Translations are applied
%     after rotations.
%   options may have the following fields:
%     open1d: if present, opens a 1d transfer function file (see IV3D
%       documentation)
%     render: by default, each frame is rendered to a PNG file. If you want
%       to suppress rendering for some frames, supply a logical vector of
%       length n_frames, and set it true for the frames you want to render.
%       This is useful if you need to re-start rendering.
%     size (default [1024 1024]): the size of the rendering window, and
%       therefore the number of pixels in each frame
%     lighting (default false): lighting is on if true, off if false

% Copyright 2011 by Timothy E Holy & Pei S Xu

  %% Parse arguments
  if (nargin < 6)
    options = struct;
  end
  n_frames = size(rotatesequence,1);
  options = default(options,...
    'render',true(1,n_frames),...
    'size',[1024 1024],...
    'lighting',false);

  if (size(translatesequence,1) ~= n_frames)
    error('rotation and translation sequences must be of the same size');
  end

  %% Open the file
  [fid,msg] = fopen(scriptname,'w');
  if (fid < 1)
    error(msg);
  end
  closeOnExit = onCleanup(@() fclose(fid));
  
  %% Preliminaries
  % Issue the command to load the data
  if ~isempty(uvf_filename)
    if ~exist(uvf_filename,'file')
      error([uvf_filename ' doesn''t exist on the current path']);
    end
    fprintf(fid,'open "%s"\n',uvf_filename);
  end
  
  % Issue the command to load a transfer function
  if isfield(options,'open1d')
    fprintf(fid,'open1d "%s"\n',options.open1d);
  end
  
  % Set the lighting
  lightingstr = {'false','true'};
  fprintf(fid,'lighting %s\n',lightingstr{1+options.lighting});
  
  
  %% Generate the individual frames
  rotateXYZ = {'rotateX','rotateY','rotateZ'};
  n_digits = ceil(log10(n_frames+1));
  pngformatstr = [pngbasename '%0' num2str(n_digits) 'd.png'];
  for i = 1:n_frames
    for j = 1:3
      if (rotatesequence(i,j) ~= 0)
        fprintf(fid,'%s %g\n',rotateXYZ{j},rotatesequence(i,j));
      end
    end
    if any(translatesequence(i,:) ~= 0)
      fprintf(fid,'translate %g %g %g\n',translatesequence(i,:));
    end
    if options.render(i)
      thisframename = sprintf(pngformatstr,i-1);
      fprintf(fid,'capturesingle "%s"\n',thisframename);
    end
  end
  
  % Cause the rendering window to stay open when complete
  fprintf(fid,'stayopen');
end
