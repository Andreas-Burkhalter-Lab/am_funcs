function indx = mdexplore(xproj,xorig,options)
% MDEXPLORE: explore multidimensional data interactively
% This has two modes of operation:
%   explore: plots points in 2d, and clicking on a point displays the
%     corresponding multidimensional data
%   cluster: plot points in 2d, and user can draw a polygon around the
%     points; the function returns the index of the encircled points.
%     (Click-drag to start polygon drawing, left click to continue, right
%     click to finish.)
%
% Syntax for 'explore' mode (does not block command line):
%   mdexplore(xproj,xorig)
%   mdexplore(xproj,xorig,options)
%
% Syntax for 'cluster' mode (blocks command line until polygon is finished):
%   indx = mdexplore(xproj,xorig)
%   indx = mdexplore(xproj,xorig,options)
%
% where
%   xproj is a 2-by-N matrix containing projections into two dimensions;
%   xorig is a cell or structure array containing the original multidimensional data;
%   options is a structure with the following fields:
%     plotfunc: a function handle which knows how to plot the
%       multidimensional data (default: @plot).
%     markersize: size of points to plot. Default depends on the mode,
%       larger points (size 9) for explore mode and smaller (size 4) for
%       cluster mode.
%     hax: specifcy some axes to use (if provided, supercedes fignum and
%       subplotnum below)
%     fignum: which figure to use (default: create a new figure)
%     subplotnum: if should be put in a subplot, what numbers to use
%     axisinfo: allows you to set axis limits to match another figure
%     axistight: if true, executes "axis tight" after calling your plotting
%       function (default false)
%     poptions: anything you want to be passed on to your plotting function
%
% An example: we will split an image into blocks, calculate the average
% red, green, blue value in each block, and display this as a scatter plot.
% When you click on a single point, this displays the corresponding block
% in another figure window (make sure you do not have overlapping windows,
% or you won't see the result of clicking on a point):
%
% % Load and display the complete image
% im = imread('peppers.png');
% figure; image(im)
% % Break into blocks
% sz = size(im);
% blocksz = floor(sz(1:2)/6);
% break1 = 0:blocksz(1):sz(1);
% break2 = 0:blocksz(2):sz(2);
% for i1 = 1:length(break1)-1
%   for i2 = 1:length(break2)-1
%     imsnip{i1,i2} = im(break1(i1)+1:break1(i1+1),break2(i2)+1:break2(i2+1),:);
%   end
% end
% imsnip = imsnip(:)';
% % Calculate mean RGB values
% rgb = zeros(3,length(imsnip));
% for i = 1:length(imsnip)
%   tmp = mean(imsnip{i},1);
%   tmp = mean(tmp,2);
%   rgb(:,i) = squeeze(tmp);
% end
% % Now call mdexplore
% mdexplore(rgb(1:2,:),imsnip,struct('plotfunc',@image))
%
% If you need to display more than 2 dimensions for your data points,
% consider using MDEXPLORE_IM instead.
%
% See also: MDEXPLORE_IM, MDEXPLORE_EPHYS.

% History: 
%   2007-07-15 (TEH)    make axistight false by default
%   2007-02-20 (RCH)    added option of giving the xorig in a structure
%                       array instead of a cell array
%   2007-01-18 (RCH)    changed so in non-polygon mode options data is stored
%                       in each point not in the axes (allows plots with
%                       different options to be plotted on top of each
%                       other)
%   2006-12-04 (RCH)    added ability to pass options to your plotfunction
%

  
% deal with options/defaults
  if (nargin < 3)
    options = struct;
  end
  options = default(options,'plotfunc',@plot,'markerColor','k');
  options = default(options,'axistight',false); % some called plotting programs set axes themselves; if so, don't want to interfere
  options = default(options,'poptions',0); % options to be passed to called plotting programs
  options = default(options,'hax',0); % 0 means not provided...
  
  mode = 'cluster';
  if (nargout == 0)
    mode = 'explore';
  end
  if ~isfield(options,'markersize')
    switch mode
      case 'cluster'
        options.markersize = 4;
      case 'explore'
        options.markersize = 4;
    end
  end

  [d,N] = size(xproj);
  if (d ~= 2)
    error('Input projection data must be a 2-by-N array');
  end

  created_new_figure = false;
  if options.hax == 0
      if isfield(options,'fignum')
          figure(options.fignum)
      else
          figure
          created_new_figure = true;
      end
      if isfield(options,'subplotnum')
          t = options.subplotnum;
          subplot(t(1),t(2),t(3))
      end
      hax = gca;
  else
      hax = options.hax;
  end
  if created_new_figure
    set(gcf,'HandleVisibility','off');
  end

  hpoints = line(xproj([1 1],:),xproj([2 2],:),...
    'MarkerFaceColor',options.markerColor,...
    'MarkerEdgeColor',options.markerColor,...
    'Marker','o',...
    'MarkerSize',options.markersize,...
    'Parent',hax);
  if isfield(options,'axisinfo')
    axis(hax,options.axisinfo);
  end
  
  % If wanting a return index, set up polygon drawing and use uiwait to
  % block until the user draws a polygon; then return the indices of the
  % points within the polygon
  if (nargout > 0)
    set(hpoints,'HitTest','off');
    set(hax,'ButtonDownFcn',@mdexplore_polygon);
    waitfor(hax,'UserData');
    pv = get(hax,'UserData');
    indx = find(inpolygon(xproj(1,:),xproj(2,:),pv(1,:),pv(2,:)));
  else
    % Otherwise, set up a callback so that when user clicks on a point, it
    % draws the corresponding md data; maximize flexibility by storing all
    % options with each point, not in the figure as a whole.
    % Set up the callback cell array
    cbarray = cell(N,1);
    for i = 1:N
      if iscell(xorig)
        cbarray{i} = {@mdexplore_plotcb,xorig{i},options};
      elseif isstruct(xorig)
        cbarray{i} = {@mdexplore_plotcb,xorig(i),options};
      else
        error('xorig must be either a cell array or a structure array')
      end 
    end
    set(hpoints,{'ButtonDownFcn'},cbarray);
    set(hax,'HitTest','off');
  end

function mdexplore_polygon(hax,eventdata)
  [pvx,pvy] = GetSelPolygon('go','b');
  set(hax,'UserData',[pvx; pvy]);

function mdexplore_plotcb(src,eventdata,plotdata,options)
%   hax = get(src,'Parent');
  plotfunc = options.plotfunc;
  poptions = options.poptions;
  if ~isstruct(poptions)
      if poptions == 0
        plotfunc(plotdata)
      else
        error('Plot options are neither a structure nor 0! What should I do?')
      end
  else
      plotfunc(plotdata,poptions);
  end
  if options.axistight
      axis tight
  end

