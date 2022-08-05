function roio = roidraw(hax,roi,done_command)
% ROIDRAW: draw regions of interest in an axis
%
% This function supports the following actions:
%   Draw a new ROI
%   Delete a ROI
%   Move a ROI
%   Resize a ROI
%   Shift the entire collection of ROIs (to handle image drift)
%
% Syntax:
%   roi_out = roidraw(hax,roi_in)
% passes in the current state of the ROIs, and waits for the user to take
% one action.
%
%   roidraw(hax,roi_in,done_command)
% instead takes a function handle as the third input, which is executed
% upon exit.  This is useful in creating GUIs which do not block the
% command line.
%
% See also: ROISTRUCT, ROIPLOT, ROIMEASURE.
  
% Copyright 2005 by Timothy E. Holy
  
  blocking = 0;
  if (nargout > 0 && nargin < 3)
    blocking = 1;
  end
  
  if isappdata(hax,'ROImarkers')
    % Get pre-existing handles for ROIs, if present
    hroi = getappdata(hax,'ROImarkers');
  else
    % Draw the ROIs & set up actions
    hroi = roiplot(hax,roi);
    set(hroi,'ButtonDownFcn','roid_roihandler');
    set(hax,'ButtonDownFcn','roid_axhandler');
  end
  
  