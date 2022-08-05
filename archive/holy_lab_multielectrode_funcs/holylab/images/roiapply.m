function intens = roiapply(ip,rois)
% ROIAPPLY: measure ROI intensities across entire files
%
% Syntax:
%   intens = roiapply(ip,rois)
%   intens = roiapply(ip,roifile)
% where
%   ip is an input IMPHYS structure array;
%   roifile is the name of a file containing saved ROI definitions;
%            OR
%   rois is the structure saved in such a file;
% and
%   intens is a nframes-by-nroi matrix giving the intensity in each frame
%     of ip.
%
% Warning: if there was drift in your prep, for the moment you'd better
% run this on the whole movie, rather than a subset. (See fixme)
%  
% See also: ROIMEASURE.
  
  % filename syntax
  if ischar(rois)
    load(rois)
  end
  
  roiTime = rois.time;
  if isscalar(roiTime)
    roiIndex = [1 1];
    roiTime = [ip(1).stacknum ip(end).stacknum];
  else
    roiIndex = [1 1:length(roiTime) length(roiTime)];
    roiTime = [ip(1).stacknum roiTime ip(end).stacknum];
  end
  
  for i = 1:length(roiTime)-1
    rangeIndex(i,:) = [find([ip.stacknum] >= roiTime(i),1,'first') ...
                       find([ip.stacknum] <= roiTime(i+1),1,'last')];
  end
  nIntervals = size(rangeIndex,1);

  for i = 1:nIntervals
    roi1 = rois.defs_orig;
    roi1.tform = rois.tform(roiIndex(i));
    roi2 = rois.defs_orig;
    roi2.tform = rois.tform(roiIndex(i+1));
    intens{i} = roimeasure({ip(rangeIndex(i,1):rangeIndex(i,2)).image},...
                           roi1,roi2,struct('skiplast',(i < nIntervals),...
                           'progress_data',[rangeIndex(i,1) length(ip)]));
  end
  intens = cat(1,intens{:});
  