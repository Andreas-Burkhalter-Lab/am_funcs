function cluster_roi_explore(R,clust,roi_def,stimulus_labels,options)
% cluster_roi_explore: connect responses after clustering to raw images
%
% After analysis, one frequently wants to go back to "raw" images to
% evaluate the quality of the data. This function allows a user to
% display processed data and go back to the original ROI definition.
%
% Syntax:
%   cluster_roi_explore(R,clust,roi_def,stimulus_labels,options)
% where
%   R is the n_stimuli-by-n_cells matrix used to represent responsiveness
%   clust is a 1-by-n_cells matrix containing the cluster number of each
%     cell
%   roi_def is a structure array, 1 entry per cell, containing the roi
%     definition.
%   stimulus_labels is a 1-by-n_stimuli cell array, containing the name
%     of each stimulus (must match that in the .imagine header)
%   options is passed on to imagesc_clusters to control the appearance of
%     the display.
%
% Usage:
%   Open stack_roi_rg for the data set you want to explore (including the
%     ROI definitions)
%   Call the function as above.
  
% Copyright 2010 by Timothy E. Holy
  
  if (nargin < 5)
    options = struct;
  end
  figure
  [himg,~,~,sortIndex] = imagesc_clusters(R,clust,options);
  set(gca,'YTick',1:length(stimulus_labels),'YTickLabel',stimulus_labels);
  set(himg,'ButtonDownFcn',@cre_click);

  function cre_click(src,evt)
    cp = get(gca,'CurrentPoint');
    cellindx = sortIndex(round(cp(1)));
    stimindx = round(cp(1,2));
    stack_roi_rg(struct('label',roi_def(cellindx).label,'display',stimulus_labels{stimindx}));
  end
end
