function msprof_temporal_register_gui(varargin)
% msprof_temporal_register_gui: correct errors in gradient timing
%
% Syntax:
%   msprof_temporal_register_gui  % runs with default file names
%   msprof_temporal_register_gui(basename)   % uses custom file names
%   msprof_temporal_register_gui(options)    % set custom parameter values
%   msprof_temporal_register_gui(basename,options)
% Different m/z values are plotted in different axes.  Each sample is
% assigned its own line color.  Look for samples that are consistently
% advanced or delayed compared to the others; if you see any, click on
% the line and, in the dialog that pops up, adjust the delay time until
% the peaks are reasonably well-aligned. Then click "Save" and close the
% figure.
%
% This function is also the entry point for "inserting" data about the
% HPLC solvent gradient into the data set.  If you want to do this,
% create a file (in the same directory as the raw data) with the name
% "set_gradients." The syntax should be as follows:
%   gc = set_gradients(flc)
% where
%   flc is a cell array of file names
%   gc is a cell array, each entry giving the gradient supplied for the
%     corresponding files in flc
% The syntax of the gradient information is the following:
%    [t1    c1; ...
%     t2    c2; ...
%     t3    c3; ...
%     ...]
% where c_i is the %B at time t_i.  The concentration is linearly
% interpolated between times, and you need to make sure you encode this
% accurately.  For example, suppose you start at 10%B, hold there for 5
% minutes, then ramp up to 80%B over 15 minutes, then hold it for the
% next 10 minutes, then drop back down to 10%B over 2 minutes.  This
% would be supplied as
%    [0 10;
%     5 10;
%     20 80;
%     30 80;
%     32 10];
% You really should plot this information to check that you've gotten it
% right, e.g.,
%    plot(m(:,1),m(:,2))
%
% See also: msprof_nnmf, msprof_choose_peak_timecourse_gui.
  
% Copyright 2010 by Timothy E. Holy
    
  options = struct;
  basename = 'msprof_ecs';
  for i = 1:length(varargin)
    if isstruct(varargin{i})
      options = varargin{i};
    elseif ischar(varargin{i})
      basename = varargin{i};
    end
  end
  
  options = default(options,...
    'files_filename',[basename '_filelist.mat'],...
    'nnmf_results_filename',[basename '_nnmf.mat'],...
    'gridsz',[5 5],...
    'index',1:24);
  if (length(options.index) > prod(options.gridsz))
    error('Grid size is insufficient for the # of "channels" you have selected');
  end
  
  pos = get(0,'ScreenSize');
  hfig = figure('Visible','off','Position',pos);
  set(hfig,'Units','normalized');
  uicontrol('Style','pushbutton','String','Save',...
    'Units','normalized',...
    'Position',[0.9 0.1 0.05 0.05],...
    'Callback',@msptrg_save);
  set(hfig,'Visible','on')
  
  nnmfdata = load(options.nnmf_results_filename);
  n_files = length(nnmfdata.tc);
  toffset = zeros(1,n_files);
  setappdata(hfig,'toffset',toffset);

  if exist('set_gradients.m','file')
    nnmfdata.gc = set_gradients(nnmfdata.file);
    nnmfdata.xaxisstr = '%B';
    nnmfdata.xc = msptrg_g2x(nnmfdata,toffset);
  else
    warndlg('Gradient data not available');
    nnmfdata.xaxisstr = 'Time (min)';
    nnmfdata.xc = nnmfdata.tc;
  end
  setappdata(hfig,'nnmfdata',nnmfdata)
  
  col = zeros(n_files,3);
  for i = 1:n_files
    col(i,:) = unique_color(i+1,n_files+1);  % avoid using black as first color
  end
  options = default(options,'color',col);
  setappdata(hfig,'options',options);
  
  msptrg_update(hfig)
  
function msptrg_update(hfig)
  options = getappdata(hfig,'options');
  nnmfdata = getappdata(hfig,'nnmfdata');
  n_files = length(nnmfdata.tc);
  keepFlag = cell(1,n_files);
  for j = 1:n_files
    keepFlag{j} = ~isnan(nnmfdata.xc{j});
  end
  for i = 1:length(options.index)
    subplot(options.gridsz(1),options.gridsz(2),i);
    cla
    for j = 1:n_files
      line(nnmfdata.xc{j}(keepFlag{j}), ...
	   nnmfdata.nnmfresults(options.index(i)).iProf{j}(keepFlag{j}),...
	   'Color',options.color(j,:),...
	   'ButtonDownFcn',@(sender,event) msptrg_click(sender,event,j));
    end
  end
      
function msptrg_click(sender,~,index)
  hfig = get_parent_fig(sender);
  nnmfdata = getappdata(hfig,'nnmfdata');
  toffset = getappdata(hfig,'toffset');
  answer = inputdlg('Enter the desired delay (min) (use - for advance)','Enter delay',1,{num2str(toffset(index))});
  if ~isempty(answer)
    tmp = str2double(answer);
    if ~isnan(tmp)
      toffset(index) = tmp;
      setappdata(hfig,'toffset',toffset);
      nnmfdata.xc = msptrg_g2x(nnmfdata,toffset);
      setappdata(hfig,'nnmfdata',nnmfdata);
      msptrg_update(hfig)
    end
  end
  
function xc = msptrg_g2x(nnmfdata,toffset)
  n_files = length(nnmfdata.tc);
  gc = nnmfdata.gc;
  % For each file, find the beginning and end of the monotonic rise
  % This is copied from resample_gradient
  gsnipc = cell(1,n_files);
  indexsnip = zeros(2,n_files);
  for k = 1:n_files
    g = gc{k};
    dg = diff(g,1,1);
    indxstart = find(dg(:,2) > 0,1,'first');
    if isempty(indxstart)
      error('There is no gradient for this file');
    end
    n = find(dg(indxstart:end,2) < 0,1,'first');
    if isempty(n)
      indxend = size(dg,1);
    else
      indxend = n+indxstart-1;
    end
    while(indxend > indxstart && g(indxend,2) == g(indxend-1,2))
      indxend = indxend-1;
    end
    g = g(indxstart:indxend,:);
    gsnipc{k} = g;
    indexsnip(:,k) = [indxstart; indxend];
  end
  xc = cell(1,n_files);
  for k = 1:n_files
    xc{k} = interp1(gsnipc{k}(:,1),gsnipc{k}(:,2),nnmfdata.tc{k} - toffset(k));
  end
  
function msptrg_save(sender,~)
  hfig = get_parent_fig(sender);
  nnmfdata = getappdata(hfig,'nnmfdata');
  nnmfdata.randKey = rand;
  options = getappdata(hfig,'options');
  toffset = getappdata(hfig,'toffset');
  nnmfdata.toffset = toffset;
  save(options.nnmf_results_filename,'-struct','nnmfdata');
  