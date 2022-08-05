function track_tube_traversals(options)
% TRACK_TUBE_TRAVERSALS: GUI for quantifying video, entrances and exits to tubes
%
% This function is specifically designed for mouse behavioral assays,
% detecting when a mouse enters or exits one or more tubes. The user
% picks a time range of the movie to analyze, and then draws circles at
% the entrance or exit of one or more tubes (and optionally, one or more
% control areas).  Intensity is then quantified in the circles and in
% annuli centered on the circles, but only in the upper half of the
% video.
%
% Syntax:
%   track_tube_traversals
% This leads the user through all aspects: selecting files, selecting
% time ranges, drawing circles, and runs the quantification.
%
%   track_tube_traversals(struct('files',{{'MC21_base_5.27.10.avi'}},'trange',{{[1500
%   5000]}}))
% This syntax lets you specify a particular file (note you could supply a
% list of files) and time range programmatically.
%
% You can control the radius of the outer annulus using the field
% 'rfactor', which defaults to sqrt(2).
% 
% See also: movie_quantify_intensity.
  
% Copyright 2010 by Timothy E. Holy
  
  if (nargin < 1)
    options = struct;
  end
  options = default(options,'rfactor',sqrt(2),'quantfunc',@(p)prctile(p,10));
  
  %% Pick movie files
  if isfield(options,'files')
    files = options.files;
  else
    files = UIGetFiles('*','Please pick your movie files');
  end
  n_files = length(files);
  trange_byfile = cell(1,n_files);
  circparams_byfile = cell(1,n_files);
  moviesz = zeros(n_files,2);
  
  %% Select the analysis time span and draw circles
  if ~exist('ttt_allparams.mat','file')
      for fileIndex = 1:n_files
        videofile = files{fileIndex};
        % Open the movie
        mov=mmreader_ffmpeg(videofile,'calibrate',true);
        im = read(mov,15);  % see if this fixes occasional crashes

        moviesz(fileIndex,:) = [mov.Height mov.Width];
        if isfield(options,'trange')
          trange = options.trange{fileIndex};
        else
          trange = select_movie_times_gui(mov,2,'start','stop');
        end
        trange_byfile{fileIndex} = trange;

        % Read & display the initial frame
        im = readAtTime(mov,trange(1));
        hfig = figure;
        himg = imshow(im);
        % Have the user draw circles
        circparams = add_circles(himg);
        close(hfig)
        circparams_byfile{fileIndex} = circparams;
        delete(mov)
      end
      save('ttt_allparams','circparams_byfile','trange_byfile','moviesz');
      drawnow
  else
      load('ttt_allparams');
  end
  
  %% Run the quantification
  % Convert circles to pixels
  for fileIndex = 1:n_files
    circles = circparams_byfile{fileIndex};
    % Add circles that are of larger size
    n_circles = size(circles,1);
    circles = [circles; circles]; %#ok<AGROW>
    circles(n_circles+1:end,3) = circles(1:n_circles,3)*options.rfactor;
    circles = [round(circles(:,1:2)) ceil(circles(:,3))];
    blobs = cell(1,2*n_circles);
    for i = 1:n_circles
      blobs{i} = ttt_circ2pixels(circles(i,1),circles(i,2),circles(i,3),moviesz(fileIndex,:));
      tmp = ttt_circ2pixels(circles(i+n_circles,1),circles(i+n_circles,2),circles(i+n_circles,3),moviesz(fileIndex,:));
      blobs{i+n_circles} = setdiff(tmp,blobs{i});
    end
    trange = trange_byfile{fileIndex};
    % Re-open the movie
    videofile = files{fileIndex};
    mov=mmreader_ffmpeg(videofile,'calibrate',true);
    [intensity,frameTime] = movie_quantify_intensity(mov,trange,blobs,options.quantfunc); %#ok<NASGU>
    save(replace_extension(videofile,'_intensity.mat'), 'videofile','trange','circles','intensity','frameTime')
    delete(mov)
  end
  delete('ttt_allparams.mat');
end

function indx = ttt_circ2pixels(x,y,r,sz)
  xl = (-r:r)+x;
  yl = (-r:0)+y;
  [X,Y] = ndgrid(xl,yl);
  keepFlag = (X-x).^2 + (Y-y).^2 <= r^2;
  indx = sub2ind(sz,Y(keepFlag),X(keepFlag));
end
