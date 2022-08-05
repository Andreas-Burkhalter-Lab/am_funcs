function im = imphysfetch(ip, options)
% IMPHYSFETCH: return image data
% Syntax:
%   im = imphysfetch(ip, options)
% where
%   ip is an IMPHYS structure;
%   im is the image contained in this structre.
%   options: current valid fields are:
%      show_progress: 1 if show progress bar
%
% If ip is a structure array, then im is a cell array of images.
%
% This function works for both real and virtual (vimage) images.
  
% Copyright 2005 by Timothy E. Holy
  
  if (length(ip) == 0)
    im = [];
    return;
  end
  
  if(nargin==1)
     options=[];
  end
  
  im = {ip.image};
  for i = 1:length(im)
    if isa(im{i},'vimage')
      im{i} = eval(im{i});    % We have to evaluate the image
    end
    
    if(isfield(options, 'show_progress') && options.show_progress==1)
      progress_bar(struct('progress', i, 'max', length(im), 'what', ['loading ' num2str(length(im)) ' frames ...']));
    end
  end
  if (length(im) == 1)
    im = im{1};
  end