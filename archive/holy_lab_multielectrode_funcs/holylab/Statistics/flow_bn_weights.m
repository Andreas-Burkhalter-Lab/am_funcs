function [yfinal,l2final,w,ystart,map,options] = flow_bn_weights(x,lgroups,varargin)
% FLOW_BN_WEIGHTS: flow points by balanced neighborhood mean shift, and report weights associated with each template
%
% Syntax:
% Basic usage: runs as many probe points as it deems necessary
%   [yfinal,l2final,w,ystart,options] = flow_bn_prob(x,lgroups,options)
% To run a fixed number n of probe points:
%   [yfinal,l2final,w,ystart,map,options] = flow_bn_prob(x,lgroups,n,options)
% To run an additional number n of probe points, and consolidate with
% previous results:
%   [yfinal,l2final,w,options] = flow_bn_prob(x,lgroups,n,yfinal,l2final,w,options)
%
% where
%   x is a d-by-N matrix of N data points in d dimensions (each column is a
%     single data point)
%   lgroups is a vector of length d that describes how the different
%     dimensions of x are to be linked.  Each coordinate is assigned a
%     "group number", and coordinates with the same group number share a
%     common scaling in the metric.  For example, if you want the metric to
%     scale differently for each coordinate, set lgroups = 1:d (assigning
%     each coordinate to a different group). If you want it to be symmetric
%     for all coordinates, set it to lgroups = ones(1,d) (thus assigning all
%     coordinates to the same group).
%
% On output:
%   yfinal is a d-by-k matrix, giving the locations of k unique "fixed"
%     points of the flow. These can be thought of as "templates" or the
%     "most typical element" in a region.
%   l2final is a d-by-k matrix, giving the metric scalings associated with each
%     column of yfinal at the fixed point
%   w is a k-by-N matrix specifying the weight assigning each data point to
%     the different fixed points.  The weights are not normalized (they
%     approximately represent a total "amount of sampling," so sum(w)>1
%     means that the point was probably sampled fairly extensively)
%   ystart is a d-by-n_collected matrix holding the starting locations of
%     each probe point tested
%   map is a 1-by-n_collected vector indicating the index of the point in
%     yfinal to which each starting point converged.
%
% The options structure may have the following fields:
% Fields that affect iteration and convergence:
%   tmax: the maximum time, in seconds, that you're willing to devote to
%     this computation. If set, this trumps all other convergence criteria
%     (when the maximum time is exceeded).
%   n_probes_mode: 'average_coverage' or 'min_coverage'. This option
%     determines how many probe points are used. (If you explicitly specify
%     n, this option is ignored.)  'min_coverage' requires that each
%     point have a summed weight bigger than wmin_quit. 'average_coverage'
%     requires that the typical coverage (defined as 1/mean(1./w)) exceeds
%     wmin_quit.  For a given wmin_quit, 'min_coverage' is much more strict
%     than 'average_coverage'. Default value: 'average_coverage'
%   wmin_quit: the threshold on weights for convergence. Details are
%     described under 'n_probes_mode'
%   plot_coverage_progress_interval (default 1): If set to n, this plots
%     "coverage" information every n iterations. Setting this to inf to
%     disable plotting entirely. If your data set is fairly small,
%     plotting can be the most time-consuming operation, so you may want to
%     set this to some number higher than 1. The displayed figure also
%     includes a "stop" button allowing you to gracefully terminate
%     iteration.
%   next_probe_mode: 'random', 'lowest_weight', 'biasedrandom'. Affects how
%     the next point is chosen. Recommended (& default) setting:
%     'biasedrandom', which chooses points with a probability in proportion
%     to 1/w. This helps insure that poorly-covered points will eventually
%     be selected.
%   min_probes (default 20): the bare minimum of probe points to be tested
%     (subject to time limitations)
% Fields that affect how length scales l2 are chosen:
%   l2_0: default starting guess for l2 (default: all 1s)
%   l2_start_mode: 'fixed' (always uses l2_0), 'closest_template',
%     'template_with_weight'. The last two use previous experience to set the
%     initial value of l2. The latter is recommended: this previous
%     experience is used only if the new starting point received a minimum
%     weighting (specified by wmin_l2) from the given template.
%   wmin_l2: minimum weighting to set l2 from a previous fixed point
%     (otherwise l2_0 is used)
% Fields that affect memory consumption:
%   w_is_sparse (default false): if true, the output weight matrix will be
%     sparse.  This can be useful for large data sets with many clusters.
%   w_sparse_thresh: the threshold for the lowest weight that will be
%     retained (default 1/N).
% Plus, you can specify any of the fields of msams_converge1.
%
% See also: FLOW_BN_STATS.

% Copyright 2009 by Timothy E. Holy

  %% Parse arguments
  [d,N] = size(x);
  yfinal = zeros(d,0);
  l2final = zeros(d,0);
  w = [];  % define this later after we know the value of w_is_sparse
  ystart = zeros(d,0);
  map = zeros(1,0);
  n_probes_mode = 'average_coverage';
  n_to_collect = -1;
  if (isempty(varargin))
    options = struct;
  else
    if isstruct(varargin{1})
      options = varargin{1};
    else
      if ~isscalar(varargin{1})
        error(['The second argument must be a structure or the # of probe' ...
          ' points you want']);
      end
      n_probes_mode = 'fixed';
      n_to_collect = varargin{1};
    end
    carg = 2;
    while (length(varargin) >= carg)
      if isnumeric(varargin{carg})
        yfinal = varargin{carg};
        l2final = varargin{carg};
        w = varargin{carg}';
        carg = carg+3;
      else
        options = varargin{carg};
        carg = carg+1;
      end
    end
  end
  if strcmp(n_probes_mode,'fixed')
    options.n_probes_mode = 'fixed';  % fixed overrides all others
  else
    options = default(options,'n_probes_mode',n_probes_mode);
  end
  options = default(options,...
    'min_probes',20,...
    'wmin_quit',0.2,...
    'plot_coverage_progress_interval',1,...
    'next_probe_mode','biasedrandom',...
    'l2_start_mode','template_with_weight',...
    'wmin_l2',0.1,...
    'l2_0',ones(d,1),...
    'w_is_sparse',false,...
    'w_sparse_thresh',1/N);
  if isempty(w)
    % We'd like to define w as k-by-N, but sparse operations are much too slow
    % that way. So work with its transpose.
    if options.w_is_sparse
      w = sparse(N,0);
    else
      w = zeros(N,0);
    end
  end
%   lgroups = lgroups(:);
%   assignment = [];
%   save convIn x lgroups assignment options
  
  %% The overall loop
  if (n_to_collect < 1)
    n_to_collect = options.min_probes;
  end
  n_collected = 0;
  isdone = n_collected >= n_to_collect;
  tstart = tic;
  hfig_coverage = [];
  if ~isinf(options.plot_coverage_progress_interval)
    hfig_coverage = figure;
    hax_coverage = axes('Position',[0.13 0.11 0.7 0.8]);
    uicontrol(hfig_coverage,...
      'Style','pushbutton',...
      'String','Stop',...
      'Units','normalized',...
      'Position',[0.85 0.3 0.12 0.1],...
      'Callback',@fbnw_stop);
    setappdata(hfig_coverage,'stop',false);
  end
  while (~isdone)
    %% Calculate the summed weight for each point (an indicator of "coverage")
    if ~isempty(w)
      wsum = full(sum(w,2));
      wsum(wsum < 1/N) = 1/N;   % regularize at 1/N
    else
      wsum = repmat(1/N,N,1);
    end
    %% Choose the next probe point from among the points in x
    switch options.next_probe_mode
      case 'rand'
        ptIndex = round(N*rand(1,1)+0.5);
      case 'biasedrandom'
        % Choose points with probability proportional to 1/w
        ciwsum = cumsum(1./wsum);
        ri = ciwsum(end)*rand(1,1);
        ptIndex = find(ri < ciwsum,1,'first');
        if (isempty(ptIndex))
          ptIndex = N;
        end
      case 'lowest_weight'
        [minw,ptIndex] = min(wsum);
      otherwise
        error(['next_probe_mode ' options.next_probe_mode ' not recognized']);
    end
    if (mod(n_collected,options.plot_coverage_progress_interval) == 0)
      axes(hax_coverage)
      semilogy(wsum,'x')
      line([1 N],[1 1]/mean(1./wsum),'Color','r')
      title(sprintf('Coverage progress (with %d probe points)',n_collected))
      ylabel('Summed weight')
      xlabel('Point index')
      set(gca,'XLim',[1 N],'TickDir','out')
      drawnow
    end
    y = x(:,ptIndex);
    ystart(:,end+1) = y; %#ok<AGROW>
    
    %% Choose the starting value of l2
    switch options.l2_start_mode
      case 'fixed'
        l2 = options.l2_0;
      case 'closest_template'
        if isempty(yfinal)
          l2 = options.l2_0;
        else
          dy = yfinal - repmat(y,1,size(yfinal,2));
          R2 = sum(dy.^2 ./ l2final,1);
          [minDist,minIndex] = min(R2);
          l2 = l2final(:,minIndex);
        end
      case 'template_with_weight'
        if isempty(w)
          l2 = options.l2_0;
        else
          [maxW,index] = max(w(ptIndex,:));
          if (maxW > options.wmin_l2)
            % There is support for this being a useful template
            l2 = l2final(:,index);
          else
            % We can't trust that this is a useful template, go with the
            % default
            l2 = options.l2_0;
          end
        end
      otherwise
        error(['l2_start_mode ' options.l2_start_mode ' not recognized']);
    end
    
    %% Flow the probe point
%     save -append convIn y l2
    [nbrInfo,status,minDist] = msams_converge1(x,y,l2,lgroups(:),options);
    if any(isnan(minDist) | isinf(minDist))
      keyboard
    end
    if (~status.converged)
      warning('bn:convergence','A probe point failed to converge (just choosing a new point)')
      continue
    end
    
    %% Determine whether this point has been seen before, and then update
    %the catalog of convergence locations and weights appropriately
    yf = nbrInfo.nbrhoodMean;
    thisw = exp(-minDist(:).^2/2);
    if options.w_is_sparse
      si = find(thisw > options.w_sparse_thresh);
    end
    [iy,index] = intersect(yfinal',yf','rows');
    if isempty(iy)
      % It's a new point
      yfinal(:,end+1) = yf; %#ok<AGROW>
      l2final(:,end+1) = nbrInfo.nbrhoodMSDisp; %#ok<AGROW>
      if options.w_is_sparse
        w(si,end+1) = thisw(si); %#ok<AGROW>
      else
        w(:,end+1) = thisw; %#ok<AGROW>
      end
      map(end+1) = size(yfinal,2); %#ok<AGROW>
    else
      % It's an old point
      if options.w_is_sparse
        w(si,index) = w(si,index) + thisw(si); %#ok<AGROW>
      else
        w(:,index) = w(:,index) + thisw; %#ok<AGROW>
      end
      map(end+1) = index; %#ok<AGROW>
    end
    
    %% See if we're done
    n_collected = n_collected+1;
    switch options.n_probes_mode
      case 'fixed'
        isdone = n_collected >= n_to_collect;
      case 'average_coverage'
        isdone = 1/mean(1./wsum) > options.wmin_quit;
      case 'min_coverage'
        isdone = all(wsum > options.wmin_quit);
      otherwise
        error(['n_probes_mode ' options.n_probes_mode ' not recognized'])
    end
    isdone = isdone & (n_collected >= options.min_probes);
    if isfield(options,'tmax')
      % Time trumps all
      if (toc(tstart) > options.tmax)
        isdone = true;
      end
    end
    % See if the user pressed the stop button
    if ishandle(hfig_coverage)
      isdone = isdone | getappdata(hfig_coverage,'stop');
    end
  end
  options.n_collected = n_collected;
  w = w';
  if ishandle(hfig_coverage)
    close(hfig_coverage);
  end
end

function fbnw_stop(hObject,event) %#ok<INUSD>
  hfig = get_parent_fig(hObject);
  setappdata(hfig,'stop',true);
end