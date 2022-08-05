function register_movie_apply(smm,regfilename,outfile_basename,options)
% register_movie_apply: implement registration on a movie
%
% Syntax:
%   register_movie_apply(smm)
%   register_movie_apply(smm,regfilename)
%   register_movie_apply(smm,regfilename,outfile_basename)
% where
%   smm is a stackmm object. Note that this does not have to be the same
%     file that you performed registration on; however, you will get a
%     warning if the names differ.  The main application for this is if you
%     want to optimize the deformation using a coarsened movie, and then
%     apply it to the full-size movie.  Changes in resolution are
%     automatically compensated.
%     Note: any data in the "badframes" field is ignored, because this
%     information is in the registration file.
%   regfilename is the name of the file to which registration information
%     was saved, for example by register_movie_gui,
%     register_movie_dfof_gui, or a .mat file holding the variables as
%     returned by running find_bad_frames (using the names described in
%     the help). In the latter case, you should add the variable
%     pixel_spacing to the file storing the results of find_bad_frames, to
%     make sure there are no resolution issues (e.g., if working with
%     coarsened movies).
%     If you do not supply regfilename, you will be prompted to select one.
%   outfile_basename is the base name (without extension) of the file you
%     want to write.  Both a header file (.mat format) and the .cam file
%     will be written for you.  If omitted, the default is to append
%     "_registered" after the file name.
%   options is a structure which may have the following fields:
%     rigid (default false): if true, only a global translation is applied
%       (no warping is used)
%     restrict_schedule (default []): if non-empty, the final movie is
%       "restricted" to reduce its resolution (number of pixels).
%     debug (default false): if true, a "trial" run is done, avoiding
%       writing data to disk and dumping output to the command line.
%
% See also: register_movie_gui, register_movie_dfof_gui.

% Copyright 2011 by Timothy E. Holy

  %% Argument parsing
  if (nargin < 2 || isempty(regfilename))
    [filename,pathname] = uigetfile('*.mat','Please pick the file containing registration data');
    regfilename = [pathname,filename];
  end
  if (nargin < 3 || isempty(outfile_basename))
    [~,outfile_basename] = fileparts(smm.filename);
    outfile_basename = [outfile_basename '_registered'];
  end
  if (nargin < 4)
    options = struct;
  end
  options = default(options,'debug',false,'rigid',false,'restrict_schedule',[]);
  header = smm.header;
  header.prec = 'single';  % save the data as single, so we can mark NaNs
  sz = smm.size;  
  
  %% Load the output of the registration process
  s = load(regfilename);
  %% Check for filename inconsistency
  if isfield(s,'filename')
    [~,regbase] = fileparts(char(s.filename));
    [~,smmbase] = fileparts(char(smm.filename));
    if ~strcmp(regbase,smmbase) && ~strcmp(regbase,[smmbase '_coarse'])
      answer = questdlg('Base name of smm and registration file do not agree; continue?','Yes','No');
      switch answer
        case {'','No'}
          return
      end
    end
  end
  % Extract shift & isbad, since this is stored in 2 different ways by the 2
  % GUIs
  if isfield(s,'shift')
    shift = s.shift;  % grayscale registration
  else
    shift = s.options.shift';  % dfof registration
  end 
  if isfield(s,'isbad')
    isbad = s.isbad;
  elseif isfield(s,'options') && isfield(s.options,'isbad')
    isbad = s.options.isbad;
  end  
  
  
  %% Compensate for different resolution than used in fitting
  if isfield(s,'pixel_spacing')
    scalefactor = s.pixel_spacing ./ header.pixel_spacing;
  else
    scalefactor = 1;
  end
  
  %% Prepare to apply any bad pixel data
  badpixels = smm.badpixels;
  if ~isempty(badpixels)
    badpixels = repmat(badpixels,[1 1 sz(3)]);
    fprintf('Masking out bad pixels with NaNs\n');
  end
  
  %% If we only have shift data, just use that
  % This occurs when the user only runs find_bad_frames
  if ~isfield(s,'udata')
    % Sanity check
    n_stacks = size(s.shift,2);
    if sum(header.nstacks) ~= n_stacks
      error('The number of stacks is inconsistent');
    end
    [fid,msg] = fopen([outfile_basename '.cam'],'w');
    if (fid < 0)
      error(msg);
    end
    pgoptions = struct('max',n_stacks);
    tic
    for stackIndex = 1:n_stacks
      stk = single(smm(:,:,:,stackIndex));
      if ~isempty(badpixels)
        stk(badpixels) = nan;
      end
      if any(isbad(:,stackIndex))
        bad_frames = find(isbad(:,stackIndex))';
        stk(:,:,bad_frames) = NaN;
      end
      thisshift = s.shift(:,stackIndex)';
      if ~any(isnan(thisshift))
        stk = image_shift(stk,scalefactor .* thisshift);
      else
        % Write a stack of NaNs if the shift is undefined
        stk = nan(size(stk));
      end
      if ~isempty(options.restrict_schedule)
        for rIndex = 1:size(options.restrict_schedule,1)
          stk = array_restrict(stk,logical(options.restrict_schedule(rIndex,:)));
        end
      end
      count = fwrite(fid,stk,header.prec);
      if (count < numel(stk))
        fclose(fid);
        error('Failed to write stack %d; might the disk be full?',stackIndex);
      end
      if (toc > 3 || stackIndex >= n_stacks)
        pgoptions.progress = stackIndex;
        pgoptions = progress_bar(pgoptions);
        tic
      end
    end
    fclose(fid);
    % Save the header---potentially adjust for change in image size
    header.camera = [header.camera '_registered'];
    header.width = size(stk,1);
    header.height = size(stk,2);
    header.wholeheader = '';
    header.pixel_spacing = header.pixel_spacing.*scalefactor;
    if ~isempty(options.restrict_schedule)
      header = rmfield(header,'um_per_pixel_xy');
    end
    save(outfile_basename,'header');
    return
  end
  %% Initialize stack geometry (particularly, prolongation)
  stk = smm(:,:,:,1);
  pyramid = array_restrict_schedule(sz(1:3),smm.header);
  pcoptions = register_phasecorr_initialize(stk,struct('pyramid',pyramid));
  %% Determine whether we need to scale the us to a different pixel_spacing
  if ~all(scalefactor == 1)
    fprintf('Resolutions don''t match; scaling the deformation by a factor of [%g %g %g] along each coordinate',scalefactor);
    % Scale the rigid translation
    shift = bsxfun(@times,shift,scalefactor);
    % Scale the deformation
    for stackIndex = 1:length(s.udata)
      if ~isempty(s.udata{stackIndex})
        for dimIndex = 1:3
          if ~isempty(s.udata{stackIndex}{dimIndex})
            s.udata{stackIndex}{dimIndex} = s.udata{stackIndex}{dimIndex} * scalefactor(dimIndex);
          end
        end
      end
    end
  end
  %% Prepare for temporal interpolation
  % Find stacks for which we have warp data
  haveu = cellfun(@(x) ~isempty(x),s.udata);
  udata = s.udata(haveu);
  u_stacknums = find(haveu);
  trange = [ceil(s.trange(1)) floor(s.trange(2))];
  if isfield(options,'trange')
    trange = options.trange;
  end
  thisstack = trange(1);
  % Initialize the starting u's
  uIndex = find(u_stacknums < trange(1),1,'last');
  if isempty(uIndex)
    uIndex = 0;
    u1 = [];  % u1 will contain the nearest deformation between the start of 
    %  the experiment and current stack for interpolation
  else
    u1 = register_phasecorr_prolong_fullsize(udata{uIndex},pcoptions);
  end
  if ~options.rigid && ~isempty(uIndex)
    u2 = register_phasecorr_prolong_fullsize(udata{uIndex+1},pcoptions);
    %  u2 will contain the nearest deformation between the current stack
    %  and end of the experiment
  end
  %% Open the output file
  if ~options.debug
    [fid,msg] = fopen([outfile_basename '.cam'],'w');
    if (fid < 0)
      error(msg);
    end
  end
  %% Loop over stacks and apply these warps
  pgoptions = struct('max',diff(trange));
  tic
  while (thisstack <= trange(2))
    if options.debug
      fprintf('Stack %d: ',thisstack);
    end
    if ~options.rigid
      if isempty(u1)
        u = u2;
        if options.debug
          fprintf('u = u2 (from stack %d), ', u_stacknums(uIndex+1));
        end
      elseif isempty(u2)
        u = u1;
        if options.debug
          fprintf('u = u1 (from stack %d), ', u_stacknums(uIndex));
        end
      else
        tr = u_stacknums([0 1]+uIndex);
        f = (thisstack-tr(1))/diff(tr);
        if options.debug
          fprintf('u interp %g, from stacks %d and %d, ',f,u_stacknums((0:1)+uIndex));
        else
          for dimIndex = 1:3
            u{dimIndex} = f*u2{dimIndex} + (1-f)*u1{dimIndex};
          end
        end
      end
    end
    pcoptions.shift = shift(thisstack,:);
    if ~options.debug
      stk = single(smm(:,:,:,thisstack));
      if ~isempty(badpixels)
        stk(badpixels) = nan;
      end
      if any(isbad(:,thisstack))
        bad_frames = find(isbad(:,thisstack))';
        stk(:,:,bad_frames) = NaN;
      end
      if ~any(isnan(pcoptions.shift))
        if options.rigid
          stkw = image_shift(stk,pcoptions.shift);
        else
          stkw = register_phasecorr_warp(u,stk,pcoptions);
        end
      else
        % Write a stack of NaNs if the shift is undefined
        stkw = nan(size(stk));
      end
      if ~isempty(options.restrict_schedule)
        for rIndex = 1:size(options.restrict_schedule,1)
          stk = array_restrict(stk,logical(options.restrict_schedule(rIndex,:)));
        end
      end
      count = fwrite(fid,stkw,header.prec);
      if (count < numel(stkw))
        fclose(fid);
        error('Failed to write stack %d; might the disk be full?',thisstack);
      end
    else
      fprintf('shift = [%g %g %g]\n',pcoptions.shift);
    end
    thisstack = thisstack+1;
    if (toc > 3 || thisstack >= trange(2))
      pgoptions.progress = thisstack-trange(1);
      pgoptions = progress_bar(pgoptions);
      tic
    end
    if (uIndex < length(u_stacknums) && thisstack > u_stacknums(uIndex+1))
      u1 = u2;
      uIndex = uIndex+1;
      if (uIndex < length(u_stacknums))
        pcoptions.shift = [0 0 0];
        u2 = register_phasecorr_prolong_fullsize(udata{uIndex+1},pcoptions);
      else
        u2 = [];
      end
    end
  end
  if ~options.debug
    fclose(fid);
    %% If a narrowed time span is being used, modify the output header
    if ~isequal(trange,[1 sz(4)])
        header.nstacks = diff(trange)+1;
        header.nframes = header.nstacks*header.frames_per_stack;
        keepIndex = trange(1):trange(2);
        header.stacktime = header.stacktime(keepIndex);
        header.stim_lookup = header.stim_lookup(keepIndex);
        header.trial_lookup = header.trial_lookup(keepIndex);
    end
    %% If a user has performed registration on a full resolution smm object
    % defined by multiple .imagine files, the new file must be annotated
    % as such.      
    if iscell(header.header_filename) && gt(size(header.header_filename,2),1) % stackmm denotes multiple file inputs as cell and single as char
        header.compound_cam = true;
    end   
    save(outfile_basename,'header');
  end
  