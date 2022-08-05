function make_coarse_movie(smm,restrict_schedule,basename)
% make_coarse_movie: generate a lower-resolution version of an imagine file
%
% This coarsens the stacks in a 4D (volume+time) recording, yielding a
% version that is of lower resolution but also requires less space and
% less processing time for any subsequent algorithms.  It can be useful
% to quickly get a sense of how a time-consuming analysis (e.g.,
% registration) might work out.
%
% Syntax:
%   make_coarse_movie(smm,restrict_schedule)
%   make_coarse_movie(smm,restrict_schedule,basename)
% where
%   smm is the stackmm object for the original file
%   restrict_schedule is a logical array of size n_restrictions-by-3; the
%     coarse-resolution stacks are created by successively calling
%     array_restrict in the following way:
%       for i = 1:n_restrictions
%         stk = array_restrict(stk,restrict_schedule(i,:));
%       end
%     For example,
%        restrict_schedule = [1 1 0; 1 1 0];
%     would cause a 4-fold coarsening along x and y, with no coarsening in
%     z.
%   basename is the base name (i.e., without extension) of the output
%     file. If you do not supply this input, the file will have the same
%     base name as the smm input, with '_coarse' appended to the end.
%
% A ".mat" file and a ".cam" file will be written; the .mat file is
% essentially the equivalent of the .imagine file, and the pair can be used
% to create a stackmm object.
  
% Copyright 2010-2011 by Timothy E. Holy
  
  %% Parse the inputs
  header = smm.header;
  if isfield(header,'pixel_spacing')
    pixel_spacing = header.pixel_spacing;
  elseif isfield(header,'um_per_pixel_xy')
    pixel_spacing = [[1 1]*header.um_per_pixel_xy ...
		    diff(header.piezo_start_stop)/(header.frames_per_stack-1)];
  end
  if (nargin < 3)
    [~,basename] = fileparts(smm.filename);
    basename = [basename '_coarse'];
  end
  %There will be a problem if basename is filename.  Fix me.
  sz = smm.size;
  n_restrict = size(restrict_schedule,1);
  
  %% Prepare the new header
  for i = 1:n_restrict
    pixel_spacing = pixel_spacing .* (1 + restrict_schedule(i,:));
  end
  header.pixel_spacing = pixel_spacing;
  header.camera = [header.camera '_coarsened'];
  if isfield(header,'wholeheader')
    header = rmfield(header,'wholeheader');
  end
  if isfield(header,'um_per_pixel_xy')
    header = rmfield(header,'um_per_pixel_xy');
  end
  if iscell(header.header_filename) && gt(size(header.header_filename,2),1) % stackmm denotes multiple file inputs as cell and single as char
    header.compound_cam = true;
  end
  header.prec = 'single';
  % Don't save it yet, because we want to find out the final size
  
  % Get info on bad pixels
  badpixels = smm.badpixels;
  if ~isempty(badpixels)
    badpixels = repmat(badpixels,[1 1 sz(3)]);
  end
  
  %% Write the movie
  fn = [basename '.cam'];
  [fid,msg] = fopen(fn,'w');
  if (fid < 0)
    error(msg)
  end
  fprintf('Stack: ')
  for stackIndex = 1:sz(4)
    fprintf('%d...',stackIndex);
    stk = single(smm(:,:,:,stackIndex));
    if ~isempty(badpixels)
      stk(badpixels) = nan;
    end
    for i = 1:n_restrict
      stk = array_restrict(stk,logical(restrict_schedule(i,:)));
    end
    if (stackIndex == 1)
      header.height = size(stk,2);
      header.width = size(stk,1);
      header.depth = size(stk,3);
      header.frames_per_stack = size(stk,3);
      save(basename,'header');
    end
    count = fwrite(fid,stk,header.prec);
    if (count < numel(stk))
      fclose(fid);
      error('Error writing .cam file');
    end
  end
  fclose(fid);
  fprintf('\n');
  
  