function register_run(smm,options)
% REGISTER_RUN: carry out multi-stack registration
% This function is designed to optimize and apply deformations to create
% registered movies. It is basically a wrapper for REGISTER_SMGRAD to
% handle the successive registration of long movies.
%
% Syntax:
%   register_run(smm,options)
% where
%   smm is a stackmm object (see STACKMM);
%   options is a structure which may have the following fields:
%     stacks: the list of stack numbers to process
%     base_stack: the stack number to which everything is registered
%       (default 2)
%     gdir_in: if present, the name of the directory from which to seek
%       data about pre-computed deformations;
%     gdir_out: if present, the name of the directory to which one saves
%       the optimized gs;
%     fit: if true, the deformation is optimized for each stack in
%       "stacks"; if not, linear interpolation is performed from whatever
%       data is available in gdir_in (default true);
%     warp_scale_in_microns: the size of the region over which the
%       warping is smooth (default 75);
%     scale: the scaling factor applied to determine the rate of
%       relaxation in each coordinate of the warp (default [1 1 0.1]);
%     outfile: if present, the name of the file used to save registered
%       movies;
%     covariant_fit: determines whether the g is optimized using
%       covariant deformation (default false);
%     covariant_warp: determines whether the registered movies are
%       generated using covariant deformation (default true);
%     imreduce: the amount of decimation used to reduce the size of the
%       images for the purpose of fitting g (not for producing the final
%       movie, which is always at full resolution) (default: [4 4 1]).
%
% Examples:
%
% To perform a relatively coarse registration on alterate members of a
% defined set of stacks, starting from scratch:
% register_run('2006-04-02-SPIM.imagine',struct('gdir_out','2006-04-02_g','stacks',stacknums(1:2:end),'warp_scale_in_microns',300))
%
% To refine g for these stacks at a higher resolution, filling in the
% intermediate ones as well:
% register_run('2006-04-02-SPIM.imagine',struct('gdir_out','2006-04-02_g','gdir_in','2006-04-02_g','stacks',stacknums,'warp_scale_in_microns',150))
%
% To calculate registered stacks for the first 8 stacks, without
% customizing the fit for each stack:
% register_run('2006-04-02-SPIM.imagine',struct('outfile','2006-04-02-SPIM_reg','gdir_in','2006-04-02_g','stacks',1:8,'fit',false))
%
% See also: REGISTER_SMGRAD, STACKMM.
  
% Copyright 2006 by Timothy E. Holy
  
  if (nargin < 2)
    options = struct;
  end
  
  %smm = stackmm(infile);
  sz = smm.size;
  if (length(sz) ~= 4)
    error('Stack dimensions not as expected');
  end
  h = smm.header;
  
  %[pathname,basename,ext] = fileparts(infile);

  if ~isfield(options,'fit')
    options.fit = true;
  end
  if ~isfield(options,'base_stack')
    options.base_stack = 2;
  end
  if ~isfield(options,'stacks')
    options.stacks = 1:sz(end);
  end
  if ~isfield(options,'covariant_fit')
    options.covariant_fit = false;
  end
  if ~isfield(options,'covariant_warp')
    options.covariant_warp = true;
  end
  if ~isfield(options,'imreduce')
    options.imreduce = [4 4 1];
  end
  if ~isfield(options,'sigma')
    pixelSpacing = [h.um_per_pixel_xy([1 1]) ...
		    diff(h.piezo_start_stop)/h.frames_per_stack];
    if ~isfield(options,'warp_scale_in_microns')
      options.warp_scale_in_microns = 75;
    end
    options.sigma = options.warp_scale_in_microns./pixelSpacing./options.imreduce;
  end
  if ~isfield(options,'scale')
    options.scale = [1 1 0.1];
  end
  if (isfield(options,'gdir_out') && ~isfield(options,'save_image'))
    % Save the warped (reduced) stacks only if there are not too many of
    % them (unless, of course, the user provides an override)
    options.save_image = length(options.stacks) < 40;
  end
  
% $$$   if ~isfield(options,'gdir_out')
% $$$     options.gdir_out = [pathname filesep basename '_g'];
% $$$     if exist(options.gdir_out,'dir')
% $$$       options.gdir_in = options.gdir_out;
% $$$     end
% $$$   end
  
  % Open the stack output file
  if (isfield(options,'outfile') && ~isempty(options.outfile))
    [fid,msg] = fopen(options.outfile,'w');
    if (fid < 0)
      error(['Error opening file ' options.outfile ': ' msg])
    end
  end
  if (options.fit && isfield(options,'gdir_out') && ~isempty(options.gdir_out))
    if ~exist(options.gdir_out,'dir')
      [success,msg] = mkdir(options.gdir_out);
      if ~success
        error(['Error creating directory ' options.gdir_out ': ' msg])
      end
    end
  end
  
  % Prepare the base image/stack
  im_size_raw = [h.width h.height h.depth];
  if options.fit
    fprintf('Preparing base stack (#%d)....',options.base_stack);
    im_base = single(smm(:,:,:,options.base_stack));
    im_base = imreduce(im_base,options.imreduce);
    psi1 = sqrt(im_base);
    fprintf('done\n')
  end

  % Determine whether we have starting guesses for g based on previous
  % work
  if isfield(options,'gdir_in')
    % Parse the available stacks
    gbasename = ['g_' rr_format_stacknum(options.base_stack,h.nstacks) '_'];
    g0source.file = dirbyname([options.gdir_in filesep gbasename '*.mat']);
    for fileIndex = 1:length(g0source.file)
      g0source.stack(fileIndex) = sscanf(g0source.file{fileIndex}(length(gbasename)+1:end),'%d');
    end
    if isempty(g0source.file)
      g0source.stack = [];
    end
    %g0source.stack = [options.base_stack g0source.stack];
  else
    g0source.stack = [];
    g0source.file = {};
    %g0source.stack = options.base_stack;
  end

  % Split the list of stacks to process into forward/backwards components
  stack_list = {options.stacks(options.stacks >= options.base_stack),...
		options.stacks(options.stacks < options.base_stack)};
  % do we want to sort them, or not? Presumably, not
  
  % Loop over forward/backward directions
  for directionIndex = 1:length(stack_list)
    prevStack = [];     % alternative initialization
    % Loop over individual stacks
    for stackIndex = 1:length(stack_list{directionIndex})
      thisStack = stack_list{directionIndex}(stackIndex);
      % Set up the initial guess g0
      % Find the nearest-neighbor stacks on either side
      allStacks = [g0source.stack prevStack options.base_stack];
      nn_less = max(allStacks(allStacks < thisStack));
      nn_more = min(allStacks(allStacks >= thisStack));
      closestStacks = [nn_less nn_more];
      if any(closestStacks == thisStack)
        % When one of them matches exactly, the other alphas (see below)
        % will be zero, so don't bother including them
        closestStacks = thisStack;
      end
%       [dist,closestStackIndex] = mindist(thisStack,allStacks);
%       closestStack = allStacks(closestStackIndex);
%       fprintf('For stack %d, the stack used for initialization was %d\n',...
%         thisStack,closestStack);
      fprintf('For stack %d, the closest stack(s) used for initialization are ',...
        thisStack);
      fprintf('%d ',closestStacks);
      fprintf('\n');
      % Get the g's for the nearest-neighbor stacks
      for nbrIndex = 1:length(closestStacks)
        if (closestStacks(nbrIndex) == prevStack)
          g0tmp{nbrIndex} = g;          % the closest was the previously-processed stack
        elseif (closestStacks(nbrIndex) == options.base_stack)
          g0tmp{nbrIndex} = register_g0(size(im_base)); % default identity transform
        else
          % Load from disk
          closestStackIndex = find(g0source.stack == closestStacks(nbrIndex));
          load([options.gdir_in filesep g0source.file{closestStackIndex}]);
          g0tmp{nbrIndex} = g;
        end
      end
      % Calculate the starting g0 by linear interpolation from the
      % nearest-neighbor stacks on either side
      g0 = g0tmp{1};
      if (length(g0tmp) == 2)
        alpha = abs(closestStacks - thisStack);
        alpha = alpha([2 1])/sum(alpha);
        for dimIndex = 1:length(g0)
          g0{dimIndex} = alpha(1)*g0tmp{1}{dimIndex}+alpha(2)*g0tmp{2}{dimIndex};
        end
      end
      clear g0tmp;   % to save memory, and to make sure we don't have settings left over next time
      % Load the image data
      im_current = single(smm(:,:,:,thisStack));
      if options.fit
        im_current_r = imreduce(im_current,options.imreduce);
        % Run the registration
        fitops = struct('covariant',options.covariant_fit,...
          'scale',options.scale);
        [g,psig,mu,err] = register_smgrad(psi1,sqrt(im_current_r),g0,options.sigma,1,fitops);
        % Save the results
        if (isfield(options,'gdir_out') && ~isempty(options.gdir_out))
          vars_to_save = {'g','err'};
          if options.save_image
            img = psig.^2;
            vars_to_save{end+1} = 'img';
          end
          save([options.gdir_out filesep 'g_' ...
            rr_format_stacknum(options.base_stack,h.nstacks) '_' ...
            rr_format_stacknum(thisStack,h.nstacks)],vars_to_save{:});
        end
      else
        im_current_r = [];
      end
      clear psig img
      if (isfield(options,'outfile') && ~isempty(options.outfile))
        % Calculate file position for this stack
        n_pixels = h.height*h.width*h.depth;
        stackTotalIndex = find(options.stacks == thisStack);
        fpos = 2*n_pixels*(stackTotalIndex-1);
        rr_fseek_pad(fid,fpos);
        % Do the warping
        if (numel(im_current) > numel(im_current_r))
          g_hires = register_expandg(g,size(im_current),'iminterp');
        else
          g_hires = g;
        end
        img = register_warp(im_current,g_hires,...
          struct('covariant',options.covariant_warp));
        if any(size(img) < im_size_raw)
          % Pad img with zeros to make it the original size
          coords = mat2cell(im_size_raw,1,ones(1,length(im_size_raw)));
          img(coords{:}) = 0;
        end
        count = fwrite(fid,img,'uint16');
        if (count < numel(img))
          error(['Error writing file ' options.outfile])
        end
        clear g_hires img;   % save memory!
      end
      prevStack = thisStack;
    end
  end
  
function s = rr_format_stacknum(stacknum,n_stacks)
  n_digits = floor(log10(n_stacks))+1;
  s = sprintf(['%0' num2str(n_digits) 'd'],stacknum);
  
function rr_fseek_pad(fid,fpos)
  status = fseek(fid,fpos,'bof');
  if (status < 0)
    % See if we are trying to seek beyond the end of the file
    status = fseek(fid,0,'eof');
    if (status < 0)
      error(['Seek error on output file ' options.outfile])
    end
    fcur = ftell(fid);
    z = zeros(1,1024^2,'uint8');  % one megabyte of zeros
    while (fcur < fpos)
      n_to_write = min(length(z),fpos-fcur);
      n_written = fwrite(fid,z(1:n_to_write),'uint8');
      if (n_written < n_to_write)
        error('Error padding file with zeros');
      end
      fcur = ftell(fid);
    end
  end
