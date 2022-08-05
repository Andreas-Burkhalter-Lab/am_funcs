function stimdata = imstimcalcsvd(ip,options)
% IMSTIMCALCSVD: analyze stimulus responses by SVD
% Syntax:
%   stimdata = imstimcalcsvd(ip)
%   stimdata = imstimcalcsvd(ip,options)
% where
%   ip is an imphys structure array;
%   options is a structure which may have the following fields:
%     flush (default 0): the number of the flush valve;
%     numpre (default 4): the number of frames used in computing the
%       background (immediately before the stimulus turns on);
%     numpost (default 6): the number of frames used after the stimulus
%       turns off;
%     valvenum (default all): valves to analyze;
%     mask (default all): if present, this should be a logical (0 or 1)
%       matrix of the same size as the images; each image will be
%       multiplied by the mask before computing the SVD;
%     nsvds (default 2): number of components for which to save
%       images (note all temporal data & singular values are saved);
%     savefid: if present, saves the image data to the file pointed
%       to by savefid, and sets the 'bg' and 'images' fields (see
%       below) to file pointers instead of image data. This is
%       useful for saving memory.
%     norm_filtersize (default 4): if present, smooths each component
%       spatially with a gaussian filter of size norm_filtersize ONLY for
%       computing the normalization <u^2>/<bg * u>. Set to < 0.25 if you
%       don't want any filtering.
% and
%   stimdata is a structure array, one element per valve opening, with
%     the following fields:
%       bg: the background image
%       images: the spatial component of the difference from background
%         (only the first 2 components)
%       time: a matrix of temporal responses, where each column refers to
%         a particular component;
%       sv: the singular values;
%       norm: the vector <u^2>/<bg * u>, where u is a spatial component
%         (each element of the vector corresponds to a component). If you
%         divide "images" by norm, and multiply "time" by norm, then you
%         should produce images & temporal response traces that have been
%         matched in terms of sign & amplitude.
%       stimulus: the valve numbers for each frame in this trial.
%
% See also: IMSTIMCALCDF.
  
% Copyright 2005 by Timothy E. Holy
  
  if (nargin < 2)
    options = struct;
  end
  if ~isfield(options,'flush')
    options.flush = 0;
  end
  if ~isfield(options,'numpre')
    options.numpre = 4;
  end 
  if ~isfield(options,'numpost')
    options.numpost = 6;
  end
  if ~isfield(options,'nsvds')
    options.nsvds = 2;
  end
  if ~isfield(options,'norm_filtersize')
    options.norm_filtersize = 4;
  end
 
  stim = [ip.stimulus];
  if ~isfield(options,'valvenum')
    options.valvenum = setdiff(unique(stim),options.flush);
  end
  
  indx_on = find(stim(1:end-1) == options.flush & ...
                 stim(2:end) ~= options.flush) + 1;
  indx_off = find(stim(1:end-1) ~= options.flush & ...
                  stim(2:end) == options.flush) + 1;
  if (stim(end) ~= options.flush)
    indx_off = [indx_off length(stim)+1];
  end
  vlv_on = stim(indx_on);
  indxkeep = indexainb(vlv_on,options.valvenum);
  indx_on = indx_on(indxkeep);
  indx_off = indx_off(indxkeep);
  ntrans = length(indx_on);
  rng_bg = 1:options.numpre;
  h = fspecial('gaussian',ceil(4*options.norm_filtersize),options.norm_filtersize);
  progress_data = progress(struct('max',ntrans,'progress',0,...
                          'what','Valve transitions'));
  for i = 1:ntrans
    % Get the range of frames
    rng = max(indx_on(i)-options.numpre,1):min(indx_off(i)+options.numpost-1,length(ip));
    % Fetch the frames themselves
    sz = size(ip(1).image);
    frames = zeros(sz(1),sz(2),length(rng));
    for j = 1:length(rng)
      frames(:,:,j) = imphysfetch(ip(rng(j)));
      if isfield(options,'mask')
        frames(:,:,j) = frames(:,:,j).*options.mask;
      end
    end
    % Do the SVD
    [im,sv,time,bg] = imsvd(frames,rng_bg);
    % Preserve/save the results
    if isfield(options,'savefid')
      stimdata(i).bg = ftell(options.savefid);
      fwrite(options.savefid,bg,'float32');
      for j = 1:options.nsvds
        stimdata(i).images(j) = ftell(options.savefid);
        fwrite(options.savefid,im(:,:,j),'float32');
      end
    else
      stimdata(i).bg = bg;
      stimdata(i).images = im(:,:,1:options.nsvds); % first components
    end
    stimdata(i).imsize = size(bg);   % Necessary if loading from disk
    stimdata(i).sv = sv;
    stimdata(i).time = time;
    % Calculate the normalization factor
    for j = 1:options.nsvds
      uf = imfilter(im(:,:,j),h);
      stimdata(i).norm(j) = sum(sum(uf.^2))/sum(sum(uf.*bg));
    end
    stimdata(i).stimulus = stim(rng);
    % Update progress bar
    progress_data.progress = i;
    progress_data = progress_bar(progress_data);
  end

