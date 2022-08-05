function track_mouse_in_ROI_gui(options)
% TRACK_MOUSE_IN_ROI: GUI for behavioral experiments, picking movies and regions
%
% This function is specifically designed for mouse behavioral assays,
% and handles the user-interaction in picking files to analyze, selecting
% time ranges within the movies, and letting the user specify the ROI.
%
% Syntax:
%   track_mouse_in_ROI_gui
% This leads the user through all aspects: selecting files, selecting
% time ranges, and drawing a polygon (using impoly) to specify the ROI.
%
%   track_mouse_in_ROI_gui(struct('files',{{'MC21_base_5.27.10.avi'}},'trange',{{[20 2000]}}))
% This syntax lets you specify a particular file (note you could supply a
% list of files) and time range programmatically. The time range units are
% in seconds.
%
% You can also modify the default parameters for single-mouse tracking,
% with the following fields in the second (structure) input:
%   rmouse (default 12): the "radius" of a mouse, in pixels.
%   tskip (default 0): the amount of time (in seconds) to skip between
%     analyzed frames. 0 means to analyze every frame.
% 
% See also: track_mouse_in_ROI_analyze, impoly.
  
% Copyright 2011 by Timothy E. Holy
  
  if (nargin < 1)
    options = struct;
  end
  options = default(options,'rmouse',12,'tskip',0);
  
  %% Pick movie files
  if isfield(options,'files')
    files = options.files;
  else
    files = UIGetFiles('*','Please pick your movie files');
  end
  n_files = length(files);
  
  %% Pick the output file
  % Doing this now ensures that any permissions errors do not cause the
  % user to waste work
  [filename,pathname] = uiputfile('*','Where should I save the results?');
  if isempty(filename)
    return
  end
  % Test to see whether we can write there
  save([pathname,filename],'files');
  
  
  %% Select the analysis time span and draw the polygon
  polygon_byfile = cell(1,n_files);
  moviesz = zeros(n_files,2);
  for fileIndex = 1:n_files
    videofile = files{fileIndex};
    % Open the movie
    mov=mmreader_ffmpeg(videofile);
    mov.FrameRate = mov.calibrateFrameRate;  % fix any problems with the frame rate
  
    moviesz(fileIndex,:) = [mov.Height mov.Width];
    if isfield(options,'trange') && length(options.trange) >= fileIndex
      trange = options.trange{fileIndex};
    else
      trange = select_movie_times_gui(mov,2,'start','stop');
    end
    trange_byfile{fileIndex} = trange;
    
    % Read & display the initial frame
    im = readAtTime(mov,trange(1));
    hfig = figure('Name','Please draw ROI');
    imshow(im);
    % Have the user draw ROI
    hpoly = impoly;
    pos = getPosition(hpoly);
    polygon_byfile{fileIndex} = pos;

    % Save the results (do this incrementally in case of crashes)
    save([pathname,filename],'files','trange_byfile','polygon_byfile','moviesz','options');
    close(hfig)
    delete(mov)
  end
  drawnow
  
 