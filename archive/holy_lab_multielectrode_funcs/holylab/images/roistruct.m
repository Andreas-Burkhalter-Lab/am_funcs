function rois = roistruct
% ROISTRUCT: a structure for holding ROI information
%
% An ROISTRUCT has the following fields, each of which is a vector with
% one entry per ROI:
%   type: 'c' = circle, 's' = sphere
%   label: a unique integer
%   x: x-coordinate
%   y: y-coordinate
%   z: z-coordinate
%   xyradius: used if type is circle or sphere
%   zradius: used if type is sphere
%
% Additionally, there are the following fields, which apply to all the
% ROIs in the structure:
%   nextlabel: an integer to use for the next newly-defined ROI
%   tform: a transformation to apply to the coordinates before actually
%     using them (see MAKETFORM and CP2TFORM).
%
% See also: ROIPLOT, ROIMEASURE, ROIDRAW, MAKETFORM, CP2TFORM.
  
  rois = struct('type','','label',[],'x',[],'y',[],'z','xyradius',[],...
                'zradius',[],'nextlabel',1,'tform',[]);
  