function funch = register_gui_utilities
% REGISTER_GUI_UTILITIES: a collection of function handles often used in
% registration, particularly by the GUI programs. 
%
% Syntax:
%   funch = register_gui_utilities
% funch is a structure whose fields are function handles. You need to
% examine the code of this function to learn what each does.
%
% See also: register_movie_gui, register_movie_dfof_gui.

% Copyright 2011 by Timothy E. Holy

  funch = struct('match_u',@match_u,...
    'upc2mg',@upc2mg,...
    'umg2pc',@umg2pc,...
    'u_initial_guess',@u_initial_guess,...
    'calculate_dfof',@calculate_dfof,...
    'fetch_dfof',@fetch_dfof,...
    'choose_good_onsets',@choose_good_onsets,...
    'choose_good_firststacks',@choose_good_firststacks,...
    'plot_udiff',@plot_udiff,...
    'set_warpops',@set_warpops,...
    'get_trange',@get_trange,...
    'create_mask',@create_mask,...
    'register_rigid',@register_rigid,...
    'play_stack',@play_stack,...
    'play_movie',@play_movie,...
    'calculate_mismatch',@calculate_mismatch,...
    'plot_udiff_gui',@plot_udiff_gui,...
    'get_prototype_stacknum',@get_prototype_stacknum);

%% Exported functions---handles-independent
% These can be used as standalone functions, no GUI required.
function u = match_u(u,arg,pcoptions)
  % Two syntaxes:
  %   u = match_u(u,desired_level,pcoptions)
  % Use this to bring u to a particular size in the hierarchy of pyramids
  %
  %   u = match_u(u,u1,pcoptions)
  % Use this to ensure that u is at least as large as u1. You can then do
  %   u1 = match_u(u1,u,pcoptions)
  % and be sure that both will be of the same size.
  %
  % u and pcoptions correspond to the phasecorr suite.
  sz = cat(1,pcoptions.pyramid.sz);
  if isempty(u)
    usz = [0 0 0];
  else
    usz = size(u{1});
    usz(end+1:3) = 1;
  end
  % Determine the calling mode & set the desired level
  if iscell(arg)
    % Matching two u's to a common size
    u1 = arg;
    u1sz = size(u1{1});
    u1sz(end+1:3) = 1;
    uszmax = max([usz;u1sz]);
    eq = sz == repmat(uszmax,size(sz,1),1);
    desired_level = find(all(eq,2));
    if isempty(desired_level)
      error('u sizes do not match any level in the pyramid hierarchy');
    end
  else
    % Bringing u to a specific pyramid level
    desired_level = arg;
  end
  if isempty(u) || isempty(u{1})
    u = cell(1,3);
    for dimIndex = 1:3
      u{dimIndex} = zeros(pcoptions.pyramid(desired_level).sz,pcoptions.class);
    end
  else
    % Match to a pyramid level
    if all(usz == 1)
      for dimIndex = 1:3
        u{dimIndex} = repmat(u{dimIndex},pcoptions.pyramid(desired_level).sz);
      end
      return
    else
      eq = sz == repmat(usz,size(sz,1),1);
      current_level = find(all(eq,2));
      if isempty(current_level)
        error('Didn''t match the level');
      end
    end
    while (current_level > desired_level)
      % Prolong u
      current_level = current_level-1;
      for dimIndex = 1:3
        u{dimIndex} = array_prolong(u{dimIndex},pcoptions.pyramid(current_level).sz);
      end
    end
    while (current_level < desired_level)
      % Restrict u
      current_level = current_level+1;
      for dimIndex = 1:3
        u{dimIndex} = array_restrict(u{dimIndex},pcoptions.pyramid(current_level).restrict);
      end
    end
  end

function umg = upc2mg(u,imsz)
  % Convert phasecorr u representation to multigrid
  if iscell(u)
    n_dims = length(u);
    umg = cat(n_dims+1,u{:});
    szu = size(u{1});
    umg = bsxfun(@rdivide,umg,reshape(imsz./szu,[ones(1,n_dims) n_dims]));
  else
    umg = u;  % handle the empty case
  end

function upc = umg2pc(umg,imsz)
  % Convert multigrid u representation to phasecorr
  szu = size(umg);
  umg = bsxfun(@times,umg,reshape(imsz./szu(1:3),[1 1 1 3]));
  upc = mat2cell(umg,szu(1),szu(2),szu(3),[1 1 1]);
  upc = upc(:)';

function u = u_initial_guess(udata,thisstacknum,base_stacknum,pcoptions,warpops)
  have_udata = ~cellfun(@isempty,udata);
  have_udata(base_stacknum) = true;
  stacknums_udata = find(have_udata);
  switch warpops.initialguess
    case 'Unwarped'
      u = [];
    case 'Previous solution'
      u = udata{thisstacknum};
      if isfield(warpops,'gap_data')
        u = match_u(u,warpops.gap_data+1,pcoptions);  % bring to desired grid size
      end
    case 'Temporal interpolation'
      prev_stacknum = max(stacknums_udata(stacknums_udata < thisstacknum));
      next_stacknum = min(stacknums_udata(stacknums_udata >= thisstacknum));
      if isempty(prev_stacknum)
        u = udata{next_stacknum};
      elseif isempty(next_stacknum)
        u = udata{prev_stacknum};
      else
        u_prev = udata{prev_stacknum};
        u_next = udata{next_stacknum};
        u_prev = match_u(u_prev,u_next,pcoptions);
        u_next = match_u(u_next,u_prev,pcoptions);
        tspan = next_stacknum - prev_stacknum;
        u = cell(1,3);
        for dimIndex = 1:3
          u{dimIndex} = (thisstacknum-prev_stacknum)/tspan * u_next{dimIndex} + ...
            (next_stacknum-thisstacknum)/tspan * u_prev{dimIndex};
        end
      end
      if isfield(warpops,'gap_data')
        u = match_u(u,warpops.gap_data+1,pcoptions);
      end
  end

function dfof = calculate_dfof(smm,params,base)
  % params must have fields "prestim" and "peristim" defining the offsets
  % for background and foreground stacks, respectively
  params = default(params,'dfof_clip',[-Inf Inf],'dfof_smooth',[]);
  pre = base+params.prestim;
  keep = ~any(params.isbad(:,pre),1);
  pre = pre(keep);
  peri = base+params.peristim;
  keep = ~any(params.isbad(:,peri),1);
  peri = peri(keep);
  
  bg = mean(single(smm(:,:,:,pre)),4);
  fg = mean(single(smm(:,:,:,peri)),4);
  % Smooth if necessary
  if ~isempty(params.dfof_smooth)
    if isa(params.dfof_smooth,'function_handle')
      bg = params.dfof_smooth(bg);
      fg = params.dfof_smooth(fg);
    else
      bg = imfilter_gaussian_mex(bg,params.dfof_smooth);
      fg = imfilter_gaussian_mex(fg,params.dfof_smooth);
    end
  end
  % Calculate deltaF/F
  dfof = fg./bg - 1;
  % Clip deltaF/F if it exceeds the valid range (this helps limit the
  % impact of artifacts on registration)
  dfof(dfof < params.dfof_clip(1)) = params.dfof_clip(1);
  dfof(dfof > params.dfof_clip(2)) = params.dfof_clip(2);

function dfof = fetch_dfof(smm,params,base)
  % requires dfof_dir to be set
  filename = [params.dfof_dir filesep sprintf('%05d.mat',base)];
  if exist(filename,'file')
    load(filename);
  else
    dfof = calculate_dfof(smm,params,base);
    save(filename,'dfof')
  end
  
function [prototype,onset,ustim,trange] = choose_good_onsets(smm,options,offsetPre,offsetPeri)
  % Note: formerly this function insisted that every stack was good in the
  % entire stimulus sequence. Now, it has been changed so that it's good as
  % long as there is at least one good "pre" stack and at least one good
  % "peri" stack.  TEH 2011-09-06
  sz = smm.size;
  header = smm.header;
  trange = [1+max(0,-min(offsetPre)) sz(4)+min(0,-max(offsetPeri))];
  if isfield(header,'stim_lookup') && ~isempty(header.stim_lookup)
    [onset,ustim] = find_stimulus_start(header.stim_lookup,[],trange);
  else
    ustim = [];
    if ~isfield(options,'stack_skip')
      error('When you do not have stimulus information, you must supply the stack_skip to options.');
    end
    onset = {(trange(1):options.stack_skip:trange(2))'};
  end
  n_stimuli = length(onset);
  prototype = nan(1,n_stimuli);
  for stimIndex = 1:n_stimuli
    [~,sortOrder] = sort(abs(onset{stimIndex} - options.base_stacknum));
    stacknum = onset{stimIndex}(sortOrder);
    isfirst = true;
    for trialIndex = 1:length(stacknum)
      isbadPre = any(options.isbad(:,stacknum(trialIndex)+offsetPre),1);
      isbadPeri = any(options.isbad(:,stacknum(trialIndex)+offsetPeri),1);
      if all(isbadPre) || all(isbadPeri)
        stacknum(trialIndex) = NaN;
      else
        if isfirst
          % It's the first good stack, save it as the prototype
          prototype(stimIndex) = stacknum(trialIndex);
          isfirst = false;
        end
      end
    end
    if isfirst
      warndlg(['Cannot register to stimulus ' header.stim_labels{ustim(stimIndex)} ' due to bad frames']);
    end
    onset{stimIndex} = sort(stacknum(~isnan(stacknum)));   % preserve the ones that are good
  end
  
function onset_peri = choose_good_firststacks(options,stacknum,offsetPeri)
  % This function finds the first valid stack # in each valid trial
  onset_peri = zeros(size(stacknum));
  for trialIndex = 1:length(stacknum)
    stacknumPeri = stacknum(trialIndex)+offsetPeri;
    isbadPeri = options.isbad(:,stacknumPeri);
    stacknumPeri = stacknumPeri(~any(isbadPeri,1));
    onset_peri(trialIndex) = stacknumPeri(1);
  end

function [x,udiff] = calculate_udiff(udata,pcoptions)
  haveu = cellfun(@(x) ~isempty(x),udata);
  udata = udata(haveu);
%   shift = getappdata(handles.figMain,'shift');
%   shift = shift(haveu,:);
%   dshift = diff(shift,1,1);
  x = find(haveu);
  x = x(:)';
  n_u = length(x);
  udiff = zeros(3,n_u-1);
  for i = 1:n_u-1
    u1 = udata{i};
    u2 = udata{i+1};
    u1 = match_u(u1,u2,pcoptions);
    u2 = match_u(u2,u1,pcoptions);
    for dimIndex = 1:3
      udiff(dimIndex,i) = sqrt(mean((u1{dimIndex}(:) - u2{dimIndex}(:)).^2));
    end
  end

function plot_udiff(hax,udata,pcoptions,pixel_spacing)
  [x,udiff] = calculate_udiff(udata,pcoptions);
  xx = [x(1:end-1); x(2:end)];
  co = get(hax,'ColorOrder');
  for dimIndex = 1:3
    yy = repmat(udiff(dimIndex,:),2,1)*pixel_spacing(dimIndex);
    line(xx,yy,'Color',co(dimIndex,:),'Parent',hax);
  end
  
  
%% Exported functions---handles-dependent
% These must be used in a GUI that has set up appdata with the correct
% names
function ok = set_warpops(handles)
  % Set the options via the dialog
  warpops = getappdata(handles.figMain,'warpops');
  output = register_options_dialog(warpops);
  ok = true;
  if isequal(output,0)
    ok = false;
    return
  end
  % Merge the outputs into warpops
  names = fieldnames(output);  % names
  for i = 1:length(names)
    warpops.(names{i}) = output.(names{i});
  end
  setappdata(handles.figMain,'warpops',warpops);

function trange = get_trange(handles)
  % Get start/stop values
  tmp = get(handles.dragStart,'XData');
  trange = tmp(1);
  tmp = get(handles.dragStop,'XData');
  trange(2) = tmp(1);

function create_mask(handles)
  % Collect information about maximum shift across selected time range, so
  % we can make a mask that excludes the pixels that will go over the edge
  sz = getappdata(handles.figMain,'sz');
  shift = getappdata(handles.figMain,'shift');
  nanflag = isnan(shift(:,1));
  keepIndex = find(~nanflag);
  trange = get_trange(handles);
  keepFlag = keepIndex >= trange(1) & keepIndex <= trange(2);
  keepIndex = keepIndex(keepFlag);
  shift = shift(keepIndex,:);
  % Get mask parameters from user
  maskops = getappdata(handles.figMain,'maskops');
  fn = fieldnames(maskops);
  defans = struct2cell(maskops);
  for i = 1:length(defans)
    defans{i} = num2str(defans{i});
  end
  answer = inputdlg(fieldnames(maskops),'Enter mask-creation options',1,defans);
  if isempty(answer)
    return
  end
  for i = 1:length(answer)
    answer{i} = str2num(answer{i});
  end
  maskops = cell2struct(answer,fn,1);
  setappdata(handles.figMain,'maskops',maskops);
  % Create the mask
  mask = register_shift2mask(sz,shift,maskops);
  msgbox(sprintf('Masking out %d%% of the pixels',round(100*sum(mask(:) == 0)/numel(mask))));
  % Create rmg_params
  smm = getappdata(handles.figMain,'smm');
  header = smm.header;
  base_stacknum = getappdata(handles.figMain,'base_stacknum');
  base_stack = smm(:,:,:,base_stacknum);
  rmg_params = register_multigrid_options(double(base_stack),mask,struct('pixel_spacing',header.pixel_spacing));
  setappdata(handles.figMain,'rmg_params',rmg_params);
  
function register_rigid(handles)
  stacknums = getappdata(handles.figMain,'stacknums');
  selected = getappdata(handles.figMain,'selected');
  stacknums = stacknums(selected);
  % Get the base stack
  base_stacknum = getappdata(handles.figMain,'base_stacknum');
  smm = getappdata(handles.figMain,'smm');
  base_stack = single(smm(:,:,:,base_stacknum));
  % Load information about the error & shift
  err = getappdata(handles.figMain,'err');
  shift = getappdata(handles.figMain,'shift');
  % Perform the registration
  pgoptions = struct('max',length(stacknums));
  tic
  for i = 1:length(stacknums)
    thisstack = stacknums(i);
    stk = single(smm(:,:,:,thisstack));
    [stkr,thisshift] = register_rigid(base_stack,stk);
    shift(thisstack,:) = thisshift;
    imdiff = base_stack - stkr;
    err(thisstack) = nanmean(abs(imdiff(:).^2));
    if (toc > 3 || i == length(stacknums))
      pgoptions.progress = i;
      pgoptions = progress_bar(pgoptions);
      tic
    end
  end
  % Save the results
  setappdata(handles.figMain,'err',err)
  setappdata(handles.figMain,'shift',shift) 
    
function play_stack(handles)
  showDfof = false;
  if isfield(handles,'checkboxShowDfof')
    showDfof = get(handles.checkboxShowDfof,'Value');
    params = getappdata(handles.figMain,'options');
    dfof_max = str2num(get(handles.editDfofMax,'String')); %#ok<ST2NM>
    if ~isscalar(dfof_max)
      errordlg('Cannot play DeltaF/F stack until dfof_max set properly');
      return
    end
  end
  % Get the selected stack #
  stacknums = getappdata(handles.figMain,'stacknums');
  selected = getappdata(handles.figMain,'selected');
  stacknums = stacknums(selected);
  % Get the stack data
  smm = getappdata(handles.figMain,'smm');
  if showDfof
    % Selected stack
    stk = fetch_dfof(smm,params,stacknums);
    stk(stk < 0) = 0;
    stk(stk > dfof_max) = dfof_max;
    % Prototype stack
    prototype_stacknum = get_prototype_stacknum(handles,stacknums);
    base_stack = fetch_dfof(smm,params,prototype_stacknum);
    base_stack(base_stack < 0) = 0;
    base_stack(base_stack > dfof_max) = dfof_max;
  else
    % Get the base stack
    base_stacknum = getappdata(handles.figMain,'base_stacknum');
    base_stack = single(smm(:,:,:,base_stacknum));
    stk = single(smm(:,:,:,stacknums));
  end
  % Register the selected stack
  playTypes = get(handles.popupmenuPlayType,'String');
  playType = playTypes{get(handles.popupmenuPlayType,'Value')};
  switch playType
    case 'Unregistered'
      stkr = stk;
    case 'Registered, rigid'
      shift = getshift(handles,stacknums);
      if showDfof
        shift = shift - getshift(handles,prototype_stacknum);
      end
      stkr = image_shift(stk,shift);
    case 'Registered, warped'
      shift = getshift(handles,stacknums);
      if showDfof
        shift = shift - getshift(handles,prototype_stacknum);
      end
      udata = getappdata(handles.figMain,'udata');
      u = udata{stacknums};
      pcoptions = getappdata(handles.figMain,'pcoptions');
      pcoptions.shift = shift;
      stkr = register_phasecorr_warp(u,stk,pcoptions);
  end
  % Generate a movie comparing to base stack
  sz = size(stkr);
  rgb = zeros([sz(1:2) 3 sz(3)]);
  rgb(:,:,1,:) = base_stack;
  rgb(:,:,2,:) = stkr;
  rgb = double(rgb) / double(max(rgb(:)));
  m = mplay(rgb);
  set(m.hfig,'Name',sprintf('Stack %d; %s',stacknums,playType))
  
function play_movie(handles)
  frame = str2double(get(handles.editFrame,'String'));
  smm = getappdata(handles.figMain,'smm');
%   % Get the base frame
%   base_stacknum = getappdata(handles.figMain,'base_stacknum');
%   base_frame = single(smm(:,:,frame,base_stacknum));
  % Get the selected stack number(s)
  stacknums = getappdata(handles.figMain,'stacknums');
  selected = getappdata(handles.figMain,'selected');
  stacknums = stacknums(selected);
  if isempty(stacknums)
    return
  end
  regDfof = strcmp(handles.caller,'register_movie_dfof_gui');
  if regDfof
    options = getappdata(handles.figMain,'options');
    stacknums_no_prototypes = setdiff(stacknums,options.prototype);
  else
    stacknums_no_prototypes = stacknums;
  end
  % Determine what is being played
  playTypes = get(handles.popupmenuPlayType,'String');
  playType = playTypes{get(handles.popupmenuPlayType,'Value')};
  showDfof = false;
  if isfield(handles,'checkboxShowDfof')
    showDfof = get(handles.checkboxShowDfof,'Value');
    params = getappdata(handles.figMain,'options');
    dfof_max = str2num(get(handles.editDfofMax,'String')); %#ok<ST2NM>
    if ~isscalar(dfof_max)
      errordlg('Cannot play DeltaF/F movie until dfof_max set properly');
      return
    end
  end
  offset = 0;
  if regDfof && ~showDfof
    offset = 1;
  end
  % Get the shifts or deformation
  if strncmp(playType,'Registered',10)
    shift = getappdata(handles.figMain,'shift');
    shift = shift(stacknums,:);
    if any(isnan(shift(:)))
      errordlg('The shift (rigid registration) is not defined for all selected stacks');
      return
    end
    if strcmp(playType,'Registered, warped')
      udata = getappdata(handles.figMain,'udata');
      if any(cellfun(@isempty,udata(stacknums_no_prototypes)))
        errordlg('The deformation (warp registration) is not defined for all selected stacks');
        return
      end
      udata = udata(stacknums);
      pcoptions = getappdata(handles.figMain,'pcoptions');
    end
  end
  % Iteratively get stacks, register them, and select the frame
  sz = smm.size;
  if showDfof
    mov = zeros([sz(1:2) length(stacknums)],'single');
  else
    mov = zeros([sz(1:2) length(stacknums)],smm.type);
  end
  pgoptions = struct('max',length(stacknums));  % progress bar
  tic
  for i = 1:length(stacknums)
    thisstk = stacknums(i)-offset;
    if (thisstk < 1)
      continue
    end
    if showDfof
      dfof = fetch_dfof(smm,params,thisstk);
      dfof(dfof < 0) = 0;
      dfof(dfof > dfof_max) = dfof_max;
    end
    switch playType
      case 'Unregistered'
        if showDfof
          mov(:,:,i) = dfof(:,:,frame);
        else
          mov(:,:,i) = smm(:,:,frame,thisstk);
        end
      case 'Registered, rigid'
        if showDfof
          stk = dfof;
        else
          stk = smm(:,:,:,thisstk);
        end
        % Register the selected stack
        stkr = image_shift(stk,shift(i,:));
        mov(:,:,i) = stkr(:,:,frame);
        if (toc > 3 || i == length(stacknums))
          pgoptions.progress = i;
          pgoptions = progress_bar(pgoptions);
          tic
        end
      case 'Registered, warped'
        if showDfof
          stk = dfof;
        else
          stk = smm(:,:,:,thisstk);
        end
        % Register the selected stack
        u = udata{i};
        pcoptions.shift = shift(i,:);
        stkr = register_phasecorr_warp(u,single(stk),pcoptions);
        mov(:,:,i) = stkr(:,:,frame);
        if (toc > 3 || i == length(stacknums))
          pgoptions.progress = i;
          pgoptions = progress_bar(pgoptions);
          tic
        end
    end
  end
  % Play the movie
  movmax = max(mov(:));
  mov = uint8(mov*(255/double(movmax)));
  m = mplay(mov);
  set(m.hfig,'Name',sprintf('Frame %d; %s',frame,playType))
  
function prototype_stacknum = get_prototype_stacknum(handles,this_stacknum)
  smm = getappdata(handles.figMain,'smm');
  header = smm.header;
  if ~isempty(header.stim_lookup)
    this_stimulus = header.stim_lookup(this_stacknum);
  else
    this_stimulus = 1;
  end
  options = getappdata(handles.figMain,'options');
  prototype_stacknum = options.prototype(this_stimulus);
  
function plot_udiff_gui(handles)
  % Get pixel spacing so we can use physical units
  smm = getappdata(handles.figMain,'smm');
  header = smm.header;
  pixel_spacing = header.pixel_spacing;
  % Calculate and plot udiff
  udata = getappdata(handles.figMain,'udata');
  pcoptions = getappdata(handles.figMain,'pcoptions');
  cla(handles.axesUDiff);
  plot_udiff(handles.axesUDiff,udata,pcoptions,pixel_spacing);
  
function err = calculate_mismatch(handles,stacknum,mode,u)
  smm = getappdata(handles.figMain,'smm');
  options = getappdata(handles.figMain,'options');
  useDfof = ~isempty(options);
  if useDfof
    % Comparing dfof
    prototype_stacknum = get_prototype_stacknum(handles,stacknum);
    fixed = fetch_dfof(smm,options,prototype_stacknum);
    moving = fetch_dfof(smm,options,stacknum);
  else
    % Comparing raw images
    base_stacknum = getappdata(handles.figMain,'base_stacknum');
    fixed = single(smm(:,:,:,base_stacknum));
    moving = single(smm(:,:,:,stacknum));
  end
  switch mode
    case 'Unregistered'
    case 'Registered, rigid'
      shift = getshift(handles,stacknum);
      if useDfof
        shift = shift - getshift(handles,prototype_stacknum);
      end
      moving = image_shift(moving,shift);
    case 'Registered, warped'
      shift = getshift(handles,stacknum);
      if useDfof
        shift = shift - getshift(handles,prototype_stacknum);
      end
      if (nargin < 4)
        udata = getappdata(handles.figMain,'udata');
        u = udata{stacknum};
      end
      pcoptions = getappdata(handles.figMain,'pcoptions');
      pcoptions.shift = shift;
      moving = register_phasecorr_warp(u,moving,pcoptions);
    otherwise
      error('mode not recognized');
  end
  err = nanmean((fixed(:)-moving(:)).^2);
  

%% Non-exported functions
function shift = getshift(handles,stacknum)
  shift = getappdata(handles.figMain,'shift');
  shift = shift(stacknum,:);
  if any(isnan(shift))
    warndlg('Shift not defined for this stack')
  end
  shift(isnan(shift)) = 0;

  