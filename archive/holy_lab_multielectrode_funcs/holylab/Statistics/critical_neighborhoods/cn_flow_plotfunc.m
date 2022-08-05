function params = cn_flow_plotfunc(x,x0,nbrList,params)
% flow_plotfunc: display neighborhood membership in 2d during cn_flow
% Syntax:
%   params = cn_flow_plotfunc(x,x0,nbrList,params)
% where
%   x is the matrix of all data points (d-by-npts)
%   x0 is the current base point position
%   nbrList is a vector listing the columns of x that are in the
%     neighborhood
%   params is a structure with the following fields:
%     hlineNbrs (required): the line handle representing the neighbors (the
%       XData and YData will be updated with the current neighbor coords).
%       Set to empty if you do not want to plot the neighbors.
%     hlineBasePoint (required): the line handle(s) representing the
%       current base point position (the XData and YData will be set to
%       x0(1) and x0(2), respectively). Set to empty if you do not want to
%       plot the basepoint.
%     hlineBasePointTrajectory (required): the line handle(s) representing the
%       base point trajectory (the XData and YData will be extended with
%       the latest value). Set to empty if you do not want to plot the
%       basepoint trajectory. If you want to plot the x and y values as a
%       function of iteration, then define this to be a 2-vector of line
%       handles.
%     saveplot (required): if true, the plot will be saved to disk; the
%       following fields are then required:
%         hfig: the handle for the "points" (neighborhood) figure
%         imcounter: an integer the counter for the current image (will be
%           incremented by 1 before saving)
%         fmtstr: a format string, the filename will be generated with
%           sprintf(ftmstr,imcounter)
%
% Example:
%     params.hfig = figure('Color','w');  % figure for showing the neighbors
%     plot(x(1,:),x(2,:),'b.');    % plot all the points in blue
%     % Plot the neighbors in green
%     hline = line(nan,nan,'LineStyle','none','Color','g','Marker','.');
%     params.hlineNbrs = hline;
%     % Plot the current basepoint with a red 'x'
%     hline = line(nan,nan,'LineStyle','none','Color','r','Marker','x');
%     params.hlineBasePoint = hline;
%     set(gca,'Visible','off')  % don't show the axes, just show the points
%     figure    % figure for showing the basepoint trajectory
%     % Plot the trajectory in blue, the last point in red (with an x)
%     params.hlineBasePointTrajectory = plot(nan,nan,'b'); 
%     params.hlineBasePoint(2) = line(nan,nan,'LineStyle','none','Marker','x','Color','r');
%     % Set up to save each frame of the neighborhood plot
%     params.saveplot = true;
%     params.fmtstr = 'frame%04d.png'; % save as 'frame0001.png' etc
%     params.imcounter = 0;
%     % Set the updating function handle (cn_flow_by_distance, and perhaps
%     % others, looks for this)
%     params.func = @(x0,nbrList,p) cn_flow_plotfunc(x,x0,nbrList,p);
%
% See also: cn_flow_mark.

% Copyright 2011 by Timothy E. Holy

  if ~isempty(params.hlineNbrs)
    set(params.hlineNbrs,'XData',x(1,nbrList),'YData',x(2,nbrList));
  end
  if ~isempty(params.hlineBasePoint)
    set(params.hlineBasePoint,'XData',x0(1),'YData',x0(2));
  end
  if ~isempty(params.hlineBasePointTrajectory)
    if isscalar(params.hlineBasePointTrajectory)
      basepointdata = get(params.hlineBasePointTrajectory,{'XData','YData'});
      set(params.hlineBasePointTrajectory,'XData',[basepointdata{1}, x0(1)],'YData',[basepointdata{2}, x0(2)]);
    else
      basepointdata = get(params.hlineBasePointTrajectory,'YData');
      for i = 1:2
        set(params.hlineBasePointTrajectory(i),'XData',1:length(basepointdata{i})+1,'YData',[basepointdata{i}, x0(i)]);
      end
    end
  end
  drawnow
  if params.saveplot
    params.imcounter = params.imcounter+1;
    filename = sprintf(params.fmtstr,params.imcounter);
    print(params.hfig,'-dpng',filename);
  end
end