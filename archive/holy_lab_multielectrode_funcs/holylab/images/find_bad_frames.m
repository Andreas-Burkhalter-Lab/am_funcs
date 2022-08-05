function [isbad,base_stacknum,err,shift,options] = find_bad_frames(smm,options)
% find_bad_frames: compute mismatch on a per-frame basis for whole movie
%
% Syntax:
%   isbad = find_bad_frames(smm)
%   [isbad,base_stacknum,err,shift] = find_bad_frames(smm,options)
% where
%   smm is a stackmm object. If you supply badframes to start with, then
%     any additional ones are ORed with this list.
%   options is a structure which may have any of the fields:
%     base_stacknum (default mid-movie, just before a stimulus): the stack
%       to use for comparison, when computing the mismatch
%     query_base (default 'interactive'): if 'interactive', the stack is
%       displayed in a movie player and the user is prompted to
%       approve/reject it;  if 'semiauto', the stack is tested using the
%       "autonomous" badframe algorithm (see below), and if it fails the
%       user is prompted to choose a different one; if 'none', no checking
%       is performed (default 'interactive').
%     stacknums (default all): a list of stack numbers to analyze
%     coord_range (default all): a 1-by-3 cell array, each element
%       containing the range of coordinates along that dimension to "snip
%       out" of the full stackmm object.
%     update_badframes (default false):  if true, then stacknums (the set of
%       stacks to be evaluated must not be 'all' or [] and smm.badframes must
%       not be [].  This will allow for assessing the validity of each frame
%       in a badstack as well as recalculating the shift and err after
%       setting bad frames to NaNs.  One may supply a vector of stack
%       numbers with badframes from the first run to evaluate for other
%       badframes within those stacks.
%     show_progress (default true): if true, a progress bar is plotted
%     register (default true): if true, each stack is first
%       rigid-registered to the "base" stack before computing the framewise
%       mismatch
%     dx_max (default stacksize/3): the maximum displacement, in pixels,
%       permitted by rigid registration
%     registration_algorithm (default 'nancorr'):  sets the algorithm for
%       determining optimal rigid shift to apply to data. The choices are
%       'phasecorr' and 'nancorr'. Tests suggest that 'nancorr' is more
%       reliable, but slower: 'nancorr' requires more fourier transforms (5
%       to 3) on a nanpadded data set that is 2^3 times larger than phase
%       correlation, resulting in an ~16 fold increase in processing time
%     badframe_algorithm (default 'all'): choices are 'autonomous'
%       (depending on camera, looks for frames with zeros),
%       'mismatch_error' (looks for frames that don't match the base
%       stack), or 'all' (runs all algorithms and only keep frames that
%       pass all tests)
%     tfilter_size (default 5): the size of the temporal filter used to
%       determine, by median filtering, "typical" values for the error and
%       shift at any given time
%     error_factor (default 10): if the frame error exceeds the "typical"
%       error (as defined by median filtering) by this factor, it is
%       decreed a bad frame
%     shift_thresh (default [5;5;3]): if the shift differs from the
%       "typical" shift by more than this in any coordinate, the entire
%       stack is decreed bad
%     badframe_thresh (default (smm.size(1)+smm.size(2))/2)):  this threshold
%       specifies the number of pixels used in comparison of putative bad
%       frames dependent on the background bias.  Current default is tested
%       for AOB imaging with a EM Gain of 75 (backgrond bias ~900 units)
% and
%   isbad is an n_frames-by-n_stacks logical matrix, true for bad frames
%   base_stacknum is the stack number used as the base stack
%   err is an n_frames-by-n_stacks matrix, giving the mismatch error (mean
%     square pixelwise difference) for each frame in each stack. NaN means
%     the frame was missing due to registration issues
%   shift is an 3-by-n_stacks matrix giving the amount of translation
%     applied to each stack to bring it into register with the base stack.
%
% Note that isbad is calculated from err and shift, and so naturally you can
% perform a different calculation on your own (err and shift are the
% time-consuming components).
%
% See also: frame_error.

% Copyright 2011 by Timothy E. Holy

if (nargin < 2)
    options = struct;
end

header = smm.header;
sz = smm.size;
sz3 = sz(1:3);
if (sz3(3) == 1)
    sz3 = sz3(1:2);
end
n_dims = length(sz3);
if isempty(header.stim_lookup)
    prestim = 1:sum(header.nstacks);
else
    prestim = find(header.stim_lookup(1:end-1) == 0 & header.stim_lookup(2:end) > 0);
    prestim = prestim(:);
end

% If the user specified the base_stacknum directly, then by default
% disable the query
if isfield(options,'base_stacknum') && ~isempty(options.base_stacknum)
    options = default(options,'query_base','none');
end
% Set the default options
def_shift_thresh = [5;5;3];
def_shift_thresh = def_shift_thresh(1:n_dims);
options = default(options,'base_stacknum',prestim(round(length(prestim)/2)),...
    'query_base','interactive',...
    'stacknums',1:sz(4),...
    'coord_range',{[1 sz3(1)],[1 sz3(2)],[1 sz3(3)]},...
    'isbad', [],....
    'update_badframes', false,...
    'show_progress',true,...
    'register',true,...
    'registration_algorithm', 'nancorr',...
    'badframe_algorithm','all',...
    'tfilter_size',5,...
    'error_factor',10,...
    'shift_thresh',def_shift_thresh);
sz3 = cellfun(@(p) p(end)-p(1)+1,options.coord_range);
options = default(options,...
    'dx_max',sz3/3);
check_autonomous = ~isempty(strmatch(lower(options.badframe_algorithm),{'all','autonomous'},'exact'));
check_mismatch = ~isempty(strmatch(lower(options.badframe_algorithm),{'all','mismatch_error'},'exact'));
badframes = smm.badframes; % externally provided data for this data set
if ~isempty(badframes)
    if ~isempty(options.isbad)
        badframes = badframes | options.isbad;
    end
    badstacks = find(any(badframes,1));
    prestim = setdiff(prestim',find(badstacks))';
elseif ~isempty(options.isbad)
    badframes = options.isbad;
    badstacks = find(any(badframes, 1));
else
    badstacks = [];
    badframes = [];
end
if options.update_badframes
    options = default(options, 'badframe_thresh', ((sz(1)+sz(2))/2));
end
% Make the coordinate range
rng = cell(1,3);
for i = 1:3
    rng{i} = options.coord_range{i}(1):options.coord_range{i}(end);
end

% Prepare default outputs in case of early return
isbad = [];
base_stacknum = [];
err = [];
shift = [];

%% Find a "good" base stack
switch options.query_base
    case 'interactive'
        options.base_stacknum = choose_base_stack_interactively(smm,options.base_stacknum,prestim,rng);
        if isempty(options.base_stacknum)
            return
        end
    case 'semiauto'
        % Automatically check for bad frames, but prompt user if there is a
        % problem
        while true
            stk3 = smm(rng{:},options.base_stacknum);
            check_bad = autonomous_check_bad_frames(stk3,header); % for single stack not whole data set
            if ~any(check_bad)
                break
            end
            [item,ok] = listdlg('PromptString','The default base stack appears to have bad frames. Choose a different base stack #',...
                'ListString',cellstr(num2str(prestim)));
            if (ok == 0)
                return
            end
            options.base_stacknum = prestim(item);
        end
    otherwise
        fprintf('Assuming the base stack is OK.\n');
end
base_stacknum = options.base_stacknum;

badpixels = smm.badpixels;
if ~isempty(badpixels)
    badpixels = repmat(badpixels,[1 1 sz(3)]);
end

%% If registering, pre-compute fourier transform of fixed image (to save time)
fixed = double(smm(rng{:},options.base_stacknum));
if options.register
    if strcmp(options.registration_algorithm, 'nancorr')
        if ~isempty(badpixels)
            fixed(badpixels) = nan;
        end
        fixed_padded = register_translate_pad(fixed);
        fixed_fft = register_translate_nancorr(fixed_padded);
    else
        fixed_fft = fftn(fixed);
    end
end

n_stacks = length(options.stacknums);
if options.show_progress
    pgoptions = struct('max',n_stacks);
    tic
end

%% Parse each stack, looking for bad frames and (optionally) registering
isbad = false([sz3(3) n_stacks]);
err = nan([sz3(3) n_stacks]);
shift = nan(n_dims,n_stacks);
stacknums = options.stacknums;
skip_these = intersect(badstacks, stacknums);
[fid,msg] = fopen('comp.log','w');
if fid < 0
  error(msg)
end
fidClose = onCleanup(@() fclose(fid));
for i = 1:n_stacks
    stackIndex = options.stacknums(i);
    fprintf(fid,'Starting stack %d\n',stackIndex);
    drawnow('update')  % to fflush the output
    stk3 = smm(rng{:},stackIndex);
    %don't want to replace if only updating badframes such as during a second run
    if check_autonomous;
        tmp = autonomous_check_bad_frames(stk3,header);
        isbad(:,i) = tmp(:);
        clear tmp
    end
    fprintf(fid,'Past autonomous\n')
    drawnow('update')
    if options.update_badframes && check_autonomous
        [stk3_l, stk3_r] = find_adjacent_good(smm, rng, stacknums(i), skip_these);
        tmp = badstack_frame_compare(stk3, stk3_l, stk3_r, options.badframe_thresh);
        check_bad = tmp(:);
        isbad(:,i) = isbad(:,i) | check_bad(:);
        clear tmp
    end
    stk3 = double(stk3);
    if (nargout > 2 || check_mismatch)
        % Calculate the mismatch
        if options.register
            if strcmp(options.registration_algorithm,'nancorr')
                % Mark any bad pixels with NaNs
                if ~isempty(badpixels)
                    stk3(badpixels) = nan;
                end
                % Mark any bad frames marked by the user with NaNs
                if ~isempty(badframes)
                    if any(badframes(:,stackIndex))
                        stk3(:,:,badframes(:,stackIndex)) = nan;
                        isbad(:,i) = isbad(:,i) | badframes(:,stackIndex);
                    end
                end
                if any(isbad(:,i))
                    stk3(:,:,isbad(:,i)) = nan;
                end
                if ~all(isbad(:,i))
                    fprintf(fid,'About to pad\n'); drawnow('update')
                    moving_padded = register_translate_pad(stk3);
                    fprintf(fid,'About to run nancorr\n'); drawnow('update')
                    thisshift = register_translate_nancorr(fixed_fft, moving_padded,options);
                    fprintf(fid,'Done with nancorr\n'); drawnow('update')
                else
                    thisshift = [0 0 0];
                end
            elseif strcmp(options.registration_algorithm,'phasecorr')
                if ~isempty(badpixels)
                    stk3(badpixels) = 0;
                end
                % Mark any bad frames marked by the user with Zeros
                if ~isempty(badframes)
                    if any(badframes(:,stackIndex))
                        stk3(:,:,badframes(:,stackIndex)) = 0;
                        isbad(:,i) = isbad(:,i) | badframes(:,stackIndex);
                    end
                end
                if any(isbad(:,i))
                    stk3(:,:,isbad(:,i)) = 0;
                end
                if ~all(isbad(:,i))
                    fprintf(fid,'About to run register_translate\n'); drawnow('update')
                    thisshift = register_translate(fixed_fft,stk3,options);
                    fprintf(fid,'Done\n'); drawnow('update')
                else
                    thisshift = [0 0 0];
                end
                % Mark any bad pixels with NaNs for error calculation
                if ~isempty(badpixels)
                    stk3(badpixels) = nan;
                end
                % Mark any bad frames marked by the user with NaNs for error
                % calcualtion
                if ~isempty(badframes)
                    if any(badframes(:,stackIndex))
                        stk3(:,:,badframes(:,stackIndex)) = nan;
                    end
                end
                if any(isbad(:,i))
                    stk3(:,:,isbad(:,i)) = nan;
                end
            end
            fprintf(fid,'About to shift\n'); drawnow('update')
            stk3 = image_shift(stk3,thisshift);
            fprintf(fid,'Done shifting\n'); drawnow('update')
            shift(:,i) = thisshift(:);
        end
        imdiff = stk3-fixed;
        nanFlag = isnan(imdiff);
        mismatch = nansum(nansum(imdiff.^2,1),2)./sum(sum(~nanFlag,1),2);
        err(:,i) = mismatch(:);
    end
    if (toc > 3 || stackIndex == sz(end))
        pgoptions.progress = i;
        pgoptions = progress_bar(pgoptions);
        tic
    end
end

%% Check for bad frames based on mismatch
if check_mismatch
    % Look for frames that have atypically high error, or stacks that have
    % atypically large shifts. "Typical" is defined by comparison to the
    % nearby stacks, i.e., using median filtering.
    isgood = ~any(isbad,1);  % exclude any stacks already declared bad from this analysis
    errtmp = err(:,isgood);
    errfilt = medfilt1improved(errtmp',options.tfilter_size)';
    isbadtmp = errtmp > options.error_factor * errfilt;
    if options.register
        % Also look for anomalously-large shifts as a way to declare stacks as
        % being bad
        shifttmp = shift(:,isgood);
        shiftfilt = medfilt1improved(shifttmp',options.tfilter_size)';
        isbadshift = any(abs(shifttmp-shiftfilt) > repmat(options.shift_thresh(:),1,size(errtmp,2)),1);
        isbadtmp = isbadtmp | repmat(isbadshift,size(errtmp,1),1);
    end
    isbad(:,isgood) = isbadtmp;
end
end

function badframes = autonomous_check_bad_frames(stk3,header)
if strncmp(header.camera,'DV8285_BV',9)
    badframes = all(all(stk3 == 0,1),2);
else
    warning('badframes:camera','Camera not recognized, skipping bad frame check');
    badframes = false([1 1 size(stk3,3)]);
end
end

function [stk3_l, stk3_r] = find_adjacent_good(smm, rng, stacknumber, skip_these)
c = 1;
if ismember((stacknumber-1),skip_these)
    while ~isnan(c) && ismember((stacknumber-c),skip_these)
        if ~ismember((stacknumber-c-1),skip_these)
            left = stacknumber-c-1;
            c = NaN;
        else
            c = c+1;
        end
    end
else left = stacknumber - c;
end
stk3_l = smm(rng{:},left);
c = 1;
if ismember((stacknumber+1),skip_these)
    while ~isnan(c) && ismember((stacknumber+c),skip_these)
        if ~ismember((stacknumber+c+1),skip_these)
            right = stacknumber+c+1;
            c = NaN;
        else
            c = c+1;
        end
    end
else right = stacknumber + c;
end
stk3_r = smm(rng{:},right);
end

function badframes = badstack_frame_compare(stk3, stk3_l,stk3_r, thresh)
frames = size(stk3,3);
badframes = false(frames,1);
for j = 1:frames
    thisframe = double(stk3(:,:,j));
    if nansum(nansum(thisframe,2),1)==0
        badframes(j) = true;
    else
        thismin = median(min(thisframe,[],1));
        thismean = nanmean(cat(3, double(stk3_l(:,:,j)), double(stk3_r(:,:,j))),3);
        theseoff = (thisframe.^2)>(thismean+thismin).^2;
        thesetoo = ((thismean-thismin).^2)>(thisframe.^2);
        theseoff = theseoff + thesetoo;
        if sum(sum(theseoff,2),1) > thresh
            badframes(j) = true;
        end
    end
end
end
