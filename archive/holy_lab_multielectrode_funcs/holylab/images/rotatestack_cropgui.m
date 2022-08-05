function crop = rotatestack_cropgui(smm,angles,options)
% ROTATESTACK_CROPGUI: a tool for assisting choice of crop rectangles
% This function takes a series of 10 slices, rotates it, and displays the
% result. The user then drags a rectangle to choose the crop region.
%
% Syntax:
%   crop = rotatestack_cropgui(s,angles,options)
% where
%   s is a stack or a stackmm object
%   angles and options are as described in ROTATESTACK (one extra field:
%     'dimension' lets you explicitly specify the dimension around which
%     to rotate, if supplied then options.dimension = 1, 2, or 3)
% and
%   crop is a 2-by-3 matrix giving (in each column) the crop region defined
%   by the user.
%
% This routine ignores any pre-existing crop information in options.
%
% See also: ROTATESTACK.

% Copyright 2009 by Timothy E. Holy

  if isfield(options,'dimension')
    indx = options.dimension;
  else
    indx = find(angles);
    if ~isscalar(indx)
      error('There must be one nonzero angle');
    end
  end
  if isa(smm,'stackmm')
    stk = smm(:,:,:,1);
  else
    stk = smm;
  end
  ops = options;
  if isfield(ops,'crop')
    ops = rmfield(options,'crop');  % remove any pre-existing crop info
  end
  
  sz = size(stk);
  slice_coords = repmat({':'},1,3);
  n_slices = 10;
  slice_skip = ceil(sz(indx)/n_slices);
  slice_coords{indx} = 2:slice_skip:sz(indx);
  
  [sr,coords] = rotatestack(stk(slice_coords{:}),angles,ops);
  srmax = max(sr,[],indx);  % max projection of the slices
  
  nIndx = setdiff(1:3,indx);
  coordsk = coords(:,nIndx);
  
  hfig = figure;
  imshowsc(coordsk(:,2),coordsk(:,1),squeeze(srmax));
  title('Select the crop region')
  rect = getrect(hfig);
  close(hfig)
  drawnow
  
  crop = [1 1 1; sz(1:3)];
  crop(:,nIndx) = [rect(2) rect(1); rect(2)+rect(4) rect(1)+rect(3)];
end
