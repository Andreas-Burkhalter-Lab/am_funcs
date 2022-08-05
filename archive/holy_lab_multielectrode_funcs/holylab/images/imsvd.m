function [im,sv,time,bg] = imsvd(frames,bgIndex)
% IMSVD: analyze images by SVD
% Syntax:
%   [im,sv,time,bg] = imsvd(frames,bgIndex)
% where
%   frames is a 3-tensor of images, frames(:,:,k) is the kth image;
%   bgIndex is an index vector indicating the frames to use in computing
%     the background;
% and
%   im is a 3-tensor of spatial components, i.e., im(:,:,k) is the kth
%     spatial component;
%   sv is a vector of singular values;
%   time is a nframes-by-nframes matrix, where each column corresponds to
%     the temporal aspect of a component;
%   bg is the background image.
%
% See also: IMSTIMCALCSVD.
  
% Copyright 2005 by Timothy E. Holy
  
  [nrow,ncol,nframes] = size(frames);
  frames = reshape(frames,nrow*ncol,nframes);
  bg = mean(frames(:,bgIndex),2);
  frames = frames - repmat(bg,1,nframes);
  % There's a bug in the Windows version of svd, sometimes it
  % thinks it doesn't have enough memory
  % This is a workaround
  try
    [im,S,time] = svd(frames,0);
  catch
    if ~isempty(findstr(lasterr,'memory'))
      % Save data to disk, launch a new matlab process, save the
      % results, and re-load
      save frames frames
      system('matlab/r forksvd')
      load framessvd
    else
      rethrow(lasterror)
    end
  end
  sv = diag(S);
  im = reshape(im,nrow,ncol,nframes);
  bg = reshape(bg,nrow,ncol);