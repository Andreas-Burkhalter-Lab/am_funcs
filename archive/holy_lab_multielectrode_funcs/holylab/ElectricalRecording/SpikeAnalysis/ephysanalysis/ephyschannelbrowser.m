function hfig = ephysfilebrowser(ep,options)
% EPHYSFILEBROWSER: examine responses in single file
%
% Syntax:
%   ephysfilebrowser(ephys)
%   ephysfilebrowser(ephys,options)
%   hfig = ephysfilebrowser(ephys,...)  
% where
%   ephys is an EPHYS structure containing data about a single channel or
%     cell;
%   options (optional) is a structure with the following fields:
%     fieldname: 'sniptimes' or 'celltimes'. Defaults to 'celltimes' if
%       present, otherwise uses 'sniptimes'.
%     trange: a 3-vector, [pre break post], specifying the times
%       (relative to valve openings) to use in calculating DELTARATE
%       (default [-5 0 20]+0.1).
%     stimstyle: if 'both', then bars representing the mean firing rate
%       during both the pre-stimulus and peri/post-stimulus period are
%       plotted.  If 'dr', then only the difference between peri/post-
%       and pre- is shown. Default: 'both'.
%     hax: axis handles for plots (useful for GUI usage, where you want
%       to control where the data are plotted). A structure with the
%       following fields:
%         dr: handle for plotting deltarate info
%         stimulus: handle for plotting stimulus info
%         rate: handle for plotting firing rate vs. time info
%         envelope: handle for plotting envelope data
%       By default, new axes are created.  See alternative calling syntax
%         below.
%
% As an options, this function returns the handle of the created figure.
%
% See also: DELTARATE, EPHYSPLOT, EPHYSGUI.

% Copyright 2005 by Timothy E. Holy
  
  % Supply default parameters
  if (nargin < 2)
    options = struct;
  end
  if ~isfield(options,'fieldname')
    options.fieldname = 'celltimes';
    if ~isfield(ep,'sortfile')
      options.fieldname = 'sniptimes';
    end
  end
  if ~isfield(options,'trange')
    options.trange = [-5 0 20]+0.1;
  end
  if ~isfield(options,'stimstyle')
    options.stimstyle = 'both';
  end

  %
  % Input validation
  %
  if (length(ep) > 1)
    error('This works only for single files');
  end
  % Check that it's a single channel or cell
  if strcmp(options.fieldname,'sniptimes')
    if (length(ep.channels) > 1)
      error('This works only for a single channel (see ephyssubchan).');
    end
  elseif strcmp(options.fieldname,'celltimes')
    if (length(ep.cellnums) > 1)
      error('This works only for a single cell (see ephyssubcell).');
    end
  else
    error('Field type not recognized');
  end
  % Check that trange makes sense
  if (~isvector(options.trange) | length(options.trange) ~= 3)
    error('options.trange must be a 3-vector');
  end
  options.trange = options.trange(:)';
  % Check that stimstyle is specified correctly
  if isempty(strmatch(options.stimstyle,{'both','dr'},'exact'))
    error(['options.stimstyle "' options.stimstyle '" not recognized']);
  end
  
  % Determine what needs fetching
  fields = {};
  if ~isfield(ep,'stimulus')
    if isfield(ep,'stimulusfile')
      fields{end+1} = 'stimulus';
    end
  end
  if ~isfield(ep,'envelope')
    if isfield(ep,'envelopefile')
      fields{end+1} = 'envelope';
    end
  end
  if ~isfield(ep,options.fieldname)
    fields{end+1} = options.fieldname;
  end
  ep = ephysfetch(ep,fields);

  % Prepare the plot window
  if ~isfield(options,'hax')
    clf
    % Specify the relative sizes of the axes in easy-to-change format
    relsize = [1 0.4 1];   % deltarate_axis stim_axis firingrate_axis
    if isfield(ep,'envelope')
      relsize = [relsize 2];  % add the envelope_axis
    end
    relsize = relsize(end:-1:1); % Go bottom-to-top
    gap = 0.1;
    relsize = [relsize; gap+zeros(size(relsize))];
    relsize = relsize(1:end-1);
    splitpos = cumsum(relsize);
    splitpos = splitpos(1:end-1)/splitpos(end);  % Normalized positions
    hax_vect = SplitVert(splitpos);
    delete(hax_vect(2:2:end));
    hax.dr = hax_vect(1);
    hax.stimulus = hax_vect(3);
    hax.rate = hax_vect(5);
    if (length(hax_vect) > 6)
      hax.envelope = hax_vect(7);
    end
    options.hax = hax;
  end

  % Get the spike times (units: scan numbers)
  t = ep.(options.fieldname){1};
  
  % Find stimulus transition times
  stimindex = find(ep.stimulus(1,:));
  stimtime = ep.stimulus(2,stimindex);
  nstim = length(stimtime);
  
  % Calculate numbers of spikes between pairs of times surrounding a
  % stimulus
  trange_scan = round(options.trange*ep.scanrate);    % Convert to scans
  trange_scan2 = [repmat(stimtime',1,2)+repmat(trange_scan([1 2]),nstim,1); ...
                  repmat(stimtime',1,2)+repmat(trange_scan([2 3]),nstim,1)];
  npb_stim = binpairs(t,trange_scan2);
  npb_stim = [npb_stim(1:nstim);npb_stim(nstim+1:end)]';  % col1 = pre-stim,
                                                % col2 = post-stim 
  rate_stim = npb_stim ./ repmat(diff(options.trange),nstim,1);  % firing rate
  
  % Bin the firing rate into approx 1s bins
  tsplit = linspace(ep.scanrange(1),ep.scanrange(2),...
                    floor(diff(ep.scanrange)/ep.scanrate)+1);
  rate = HistSplit(t,tsplit(2:end-1))/(diff(tsplit([1 2]))/ep.scanrate);
  tcenter_sec = (tsplit(1:end-1)+tsplit(2:end))/(2*ep.scanrate);

  %
  % Do the plots
  %
  trange_sec = ep.scanrange/ep.scanrate;
  stimtime_sec = stimtime/ep.scanrate;
  % Deltarate info
  axes(options.hax.dr)
  if strcmp(options.stimstyle,'both')
    col = {'b','r'};
    for i = 1:2
      x = [stimtime_sec + options.trange(i); ...
           stimtime_sec + options.trange(i+1)];
      x = [x; x([2 1],:)];
      y = [zeros(2,nstim); repmat(rate_stim(:,i)',2,1)];
      patch(x,y,col{i});
    end
    hyl = ylabel('Avg. rates','Rotation',0,'HorizontalAlignment','right',...
                 'VerticalAlignment','middle');
  else
    dr = diff(rate_stim,1,2)';
    x = [stimtime_sec; stimtime_sec];
    y = [zeros(1,nstim); dr];
    line(x,y,'Color','k');
    hyl = ylabel('\Deltar','Rotation',0,'HorizontalAlignment','right',...
                 'VerticalAlignment','middle');
  end
  set(gca,'XLim',trange_sec);
  
  % Stimulus info
  axes(options.hax.stimulus)
  stim_off = ep.stimulus(2,stimindex+1)/ep.scanrate;
  x = [stimtime_sec; stimtime_sec; stim_off; stim_off];
  y = [zeros(1,nstim); ones(1,nstim); ones(1,nstim); zeros(1,nstim)];
  patch(x,y,'k');
  vlvstr = ep.valvelabels(ep.stimulus(1,stimindex));
  text(stimtime_sec,1.5+zeros(1,nstim),char(vlvstr));
  set(gca,'XLim',trange_sec,'Visible','off','YLim',[0 2.5]);
  hyl = ylabel('Stimulus','Rotation',0,'HorizontalAlignment','right',...
               'VerticalAlignment','middle','Visible','on');
  
  % Firing rate
  axes(options.hax.rate)
  plot(tcenter_sec,rate)
  set(gca,'XLim',trange_sec);
  hyl = ylabel('Rate','Rotation',0,'HorizontalAlignment','right',...
               'VerticalAlignment','middle');

  % Envelope
  if isfield(options.hax,'envelope')
    axes(options.hax.envelope);
    if isfield(ep,'envelope')
      ep.tag = 'all';
      pp.tags = {'all'};
      pp.fieldtoplot = 'envelope';
      pp.channelnumber = ep.channels;
      ephysplot(ep,pp);
    end
    set(gca,'XLim',trange_sec);
    hyl = ylabel('Voltage','Rotation',0,'HorizontalAlignment','right',...
                 'VerticalAlignment','middle');
  end
  
  % Clean up
  allax = [options.hax.dr options.hax.stimulus ...
           options.hax.rate];
  if isfield(options.hax,'envelope')
    allax = [allax options.hax.envelope];
  end
  set(allax(1:end-1),'XTick',[]);
  sliderwindow(allax)
  
  % Return figure handle?
  if (nargout > 0)
    hfig = gcf;
  end
  