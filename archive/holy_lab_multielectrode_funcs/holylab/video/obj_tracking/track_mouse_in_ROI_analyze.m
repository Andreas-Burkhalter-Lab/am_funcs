function track_mouse_in_ROI_analyze(filename,useroptions)
% TRACK_MOUSE_IN_ROI_ANALYZE: track a single mouse in a ROI
%
% This function is run after track_mouse_in_ROI_gui.
%
%   track_mouse_in_ROI_analyze('filename')
% where filename is the name of the file saved by track_mouse_in_ROI_gui,
% will analyze all of the movies "queued" by the GUI program.
%
%   track_mouse_in_ROI_analyze('filename',options)
% will additionally override some of the settings in the file. Any of the
% variables saved in the .mat file can be changed; for example,
%   track_mouse_in_ROI_analyze('filename',struct('tskip',1))
% will skip 1 second between analyzed frames, rather than analyzing each
% frame.
%
% The results will be saved to a MAT file with the same basename as the
% movie, with extension '.track'.
%
% You can inspect the results for veracity using
% track_mouse_in_ROI_inspect.
%
% See also: track_mouse_in_ROI_gui, track_mouse_in_ROI_inspect.

% Copyright 2011 by Timothy E. Holy

  if (nargin < 2)
    useroptions = struct;
  end
  s = load(filename);
  options = copyfields(useroptions,fieldnames(useroptions),s.options);
  n_files = length(s.files);
  for fileIndex = 1:n_files
    % Open the movie
    videofile = s.files{fileIndex};
    mov=mmreader_ffmpeg(videofile);
    mov.FrameRate = mov.calibrateFrameRate;
    % Create the mask
    sz = s.moviesz(fileIndex,:);
    p = s.polygon_byfile{fileIndex};
    mask = poly2mask(p(:,1),p(:,2),sz(1),sz(2));
    trange = s.trange_byfile{fileIndex};
    fbackground = mov.readAtTime(trange(1));
    % Create the smoothing kernel
    [X,Y] = ndgrid(1:sz(1),1:sz(2));
    XY = [X(:)'; Y(:)'];
    h = zeros(sz);
    med = sz/2; 
    keep = (X-med(1)).^2 + (Y-med(2)).^2 < options.rmouse^2;
    h(keep) = 1;
    h = h/sum(h(:));
    hfft = fftn(h);
    % Allocate storage for output
    n_frames = mov.NumberOfFrames;
    W = zeros(1,n_frames);
    center = zeros(2,n_frames);
    Cc = cell(1,n_frames);
    t = zeros(1,n_frames);
    frameIndex = 1;
    % Loop over frames
    tic;
    tstart = toc;
    while (mov.CurrentTime < trange(2))
%       fprintf('%g %g ',mov.CurrentTime,mov.CurrentTime+options.tskip);
      if options.tskip > 0
        f = mov.readAtTime(mov.CurrentTime+options.tskip);
      else
        f = mov.readNextFrame;
      end
%       fprintf('%g...',mov.CurrentTime);
      t(frameIndex) = mov.CurrentTime;
      df = fbackground(:,:,1)-f(:,:,1);
      df(~mask) = 0;
      [W(frameIndex),center(:,frameIndex),Cc{frameIndex}] = brightpatch(double(df),hfft,XY);
      if isfield(options,'display')
        display = options.display;
        if isnumeric(display)
          % Check to see if current frameIndex is one to be displayed
          if any(mov.CurrentTime == display)
            display = 'all';
          else
            display = '';
          end
        end
        if ~isempty(display)
          switch display
            case 'all'
              figure
            case 'one'
              figure(337)
            otherwise
              error('Do not recognize display option');
          end
          imshow(f)
          plot_covar_as_ellipse(gca,center([2 1],frameIndex),Cc{frameIndex}([2 1],[2 1]));
          title(num2str(mov.CurrentTime))
        end
      end
      if (toc-tstart > 5)
        fprintf('%g...',mov.CurrentTime);
        tstart = toc;
      end
      frameIndex = frameIndex+1;
    end
    fprintf('\n');
    % Chop off unused pieces
    t = t(1:frameIndex-1);
    W = W(1:frameIndex-1);
    center = center(:,1:frameIndex-1);
    Cc = Cc(1:frameIndex-1);
    % Save the data
    [pathname,basename] = fileparts(videofile);
    outname = [pathname filesep basename '.track'];
    save(outname,'t','W','center','Cc','videofile');
    % Clear the movie
    delete(mov)
  end
  
 function [W,center,C,pixIndex] = brightpatch(f,h,XY)
  if isreal(h)
    c = fftshift(ifftn(fftn(f).*conj(fftn(h))));
  else
    c = fftshift(ifftn(fftn(f).*conj(h)));
  end
  [cmax,cIndex] = max(c(:));
  % Find the connected component that includes the peak of the correlation
  keep = c > cmax/2;
  cc = bwconncomp(keep);
  cckeep = cellfun(@(s) any(s == cIndex),cc.PixelIdxList);
  if ~isempty(cckeep)
    pixIndex = cc.PixelIdxList{cckeep};
    % Compute the intensity-weighted center of mass and covariance
    w = f(pixIndex)';
    W = sum(w);
    xy = XY(:,pixIndex);
    center = sum(bsxfun(@times,xy,w),2)/W;
    dX = bsxfun(@minus,xy,center);
    C = dX*bsxfun(@times,dX,w)'/W;
  else
    W = 0;
    center = [0 0]';
    C = zeros(2,2);
    pixIndex = [];
  end
     