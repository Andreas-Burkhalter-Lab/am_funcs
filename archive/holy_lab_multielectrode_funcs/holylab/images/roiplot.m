function hroi = roiplot(hax,roi,options)
% ROIPLOT: display ROIs.
%
% Syntax:
%   roiplot(roi)
%   roiplot(hax,roi)
%   roiplot(...,options)
%   hroi = roiplot(...)
%
% where
%   hax is the axis handle used to plot the ROI (default current axis);
%   roi is an ROISTRUCT describing the ROIs
%   options is a structure which may contain the following fields:
%     linewidth: used to determine the width of lines; if this is a
%                vector, it sets each ROI's linewidth separately
%
% See also: ROISTRUCT.
  
% Copyright 2005 by Timothy E. Holy
  
  % Parse arguments
  if ~ishandle(hax)
    if (nargin > 1)
      options = roi;
    else
      options = struct;
    end
    roi = hax;
    hax = gca;
  else
    if (nargin < 3)
      options = struct;
    end
  end
  
  % todo: options = roiplot_options(options);
  
  hroi=[];
  if(isempty(roi)) 
     nroi=0;
  else
     nroi = length(roi.type);
  end
  for i = 1:nroi
    switch roi.type(i)
     case 'c'
      npts = 100;
      th = linspace(0,2*pi,npts);
      x = roi.x(i) + roi.xyradius(i)*cos(th);
      y = roi.y(i) + roi.xyradius(i)*sin(th);
      [xp,yp] = tformfwd(roi.tform,x,y);
      [textPosX, textPosY]=tformfwd(roi.tform, roi.x(i)+roi.xyradius(i), roi.y(i) );
      if(isfield(roi, 'handle') && ~isempty(roi.handle) ) % if, update position only
         if(~strcmp(roi.action{i}, 'deleted'))
            hroi(i)=roi.handle(i);
            set(hroi(i), 'XData', xp, 'YData', yp);
            set(getappdata(hroi(i), 'labelObject'), 'position', [textPosX textPosY]);
         end
      else % else, need create new line object
         % hroi(i) = line(xp,yp,'LineWidth',options.linewidth(i));
         hroi(i) = line(xp,yp, 'linewidth', 2);
         hText= text(textPosX, textPosY, num2str(roi.label(i)));
         set(hText, 'color', 'r', 'FontSize', 14); % text color is red
         
         setappdata(hroi(i), 'label', roi.label(i));
         setappdata(hroi(i), 'labelObject', hText);
      end
     otherwise
      error(['ROI type ' roi.type(i) ' not recognized'])
    end
  end

  
% todo: function  options = roiplot_options(options)
  