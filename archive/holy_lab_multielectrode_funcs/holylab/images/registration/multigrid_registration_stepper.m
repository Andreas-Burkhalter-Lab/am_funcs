function usave = multigrid_registration_stepper(infile, outfile, method, options)
% multigrid_registration_stepper multigrid_vcycle-based wrapper for various registration methods
% Syntax: usave = multigrid_registration_stepper(infile, outfiles, method, options)
% 
% Inputs: infile:  a string containing the full path to the desired input .imagine file
%                  if the infile comprises multiple .imagine files, a
%                  basefilename will be generated from the first identical
%                  fragment of the filenames (beginning with the first
%                  character and will be used in the place of a real
%                  .imagine file
%         outfile: a string containing the full path to the desired 'usave.mat' file
%         method:  a struct array with one element per "iteration" with the following fields 
%         NOTE: (some subfunctions will override supplied values)
%              **REQUIRED**
%                .func:      (no default) string identifying subfunction defined in this file (in order of application)
%                  Current accepted method.func values(2010_04_19):
%                     'fill_bad_frames' (JPM) take bad_frames and correct, attempt to identify others needing help
%                     'rigidthenmultigrid' (JPM) register_rigid, then register_multigrid_vcycle @ [0 1] relaxation  until 1% fractional error reduction
%                     'nonrigid_polish' (JPM) register_multigrid_vcycle @ [1 1] relaxation to specified stacks
%                     'linear_uinterp' (JPM) expand sparse u-values to full
%                     range using linear interpolation of u-values
%                     'write' (GFH) use previous u-values to call the
%                     warpnwrite function
%                     .cam file
%              **OPTIONAL or FUNC-SPECIFIC**
%               with method.func == 'rigidthenmultigrid'
%                .frequency: (default 'stim') string or scalar specifying the stacks to which each script should be applied
%                            can be: 'full' = all stacks
%                                    'stim' = just stacks preceding valve openings
%                                    scalar (e.g. 10) register one stack every n stacks from options.base_stack
%                .n_avg:      (default -4) (direction)(scalar) # of stacks to average at each registration point 
%
%         NOTE: each subfunction is likely to allow additional method parameters
%               Until a robust help is available, please refer to each
%               subfunction to identify additional features.
%
%         options: a struct with options to pass to register_multigrid_vcycle
%
% Outputs: usave:  struct array containing the following subfields (per method called)
%                .base_stack: scalar or vector with indices of stacks to serve as base (if vector, average will be computed)
%                .bad_stacks: vector of stack indices with known errors
%                .bad_frames: cell array of indices within each stack
%                .bad_stack_file: full path location of .cam file with stacks to replace bad_stacks
%                .stacknum:   scalar vector containing the indices of stacks to be registered
%                .u_matrix:   cell array containing (x - y - z - 3) u_data (per stack identified in stacknum field)
%                .u_err:      scalar vector same length as stacknum with error values of registered stacks
%                .orig_err:   same as u_err, but given pre-registration
%                .rmg:        rmg_params structure for this round
%                .method:     method struct used in this round
%                .multigrid_ops: options sent to register_multigrid_vcycle, if applicable
%                .lambda:     lambda value sent to register_multigrid_vcycle, if applicable 
%                .n_cycles:   (length of stacknum) # of v-cycles used in this round
%                .outfiles    cell array of full path strings with any generated output files (e.g. *.cam, usave.mat)
%
% See also register_multigrid_vcycle, register_multigrid_warp, register_multigrid_options, imcompare

% Copyright 2010 Julian P. Meeks (Tim Holy Laboratory)
% Revision history:
% 2010_04_19: Wrote it (JPM)

%% Parse infile
smm = stackmm(infile); % will throw early error w/bad input :-)

h = smm.header;
if isfield(h, 'compound_cam')
  temp = findstr('.mat', infile);
  temp = infile(1:temp-1);
  basefilename = char(temp);
elseif size(infile, 2) == 1
    temp = findstr('.imagine', infile);
    basefilename = infile(1:temp-1);
else    
  ifs = char(infile);
  count = 0;
  flag = 0;
  for i = 1:size(ifs, 2)
    if flag > 0, 
      count = count + 1;
    end
    for j = 1:(size(ifs, 1) - 1);
      if ~isequal(ifs(1,i), ifs((j+1),i))
       flag = 1;
      end
    end
      
  end
  temp = ifs(1,1:(end-count-2));
  basefilename = temp;
  end



clear h temp count flag i j ifs;
    
  
% temp = findstr('.imagine', infile);
% if ~isempty(temp)
%     basefilename = infile(1:temp-1);
% else
%     temp = findstr('.cam', infile);
%     if ~isempty(temp)
%         basefilename = infile(1:temp-1);
%     else
%         basefilename = infile;
%     end
% end


%% Parse outfile
if exist(outfile, 'file')
    load(outfile, '-mat')
    if ~exist('usave','var')
        error('supplied ''usave'' outfile is not useable by this function\n');
    elseif isfield(usave, 'infile')
        if isempty(strmatch(usave(1).infile, infile,'exact'))
            error('The supplied ''usave'' outfile is associated with a different infile.\n');
        end
    end
    for i = 1:length(usave)
        if isempty(fieldnames(usave))
            usave_length_orig = i-1;
        end
    end
    if ~exist('usave_length_orig', 'var')
        usave_length_orig = length(usave);
        usave(end+1) = fill_usave(usave(end));
    end
else
    temp = findstr('/', outfile);
    outdir = outfile(1:temp(end));
    if ~exist(outdir,'dir')
        mkdir(outdir);
    end
    usave_length_orig = 0;
end
% if necessary initialize usave variable;
if exist('usave','var')
    if usave_length_orig == 0;
        clear usave;
        usave = fill_usave(struct);
    end
else
    usave = fill_usave(struct);
end
clear temp;
%% Interpret 'method' variable
% please update this list as you add them!!!
accepted_funcs = {...
    'fill_bad_frames',...
    'rigidthenmultigrid',...
    'nonrigid_polish',... % to be written?
    'linear_uinterp',...
    'write',...
    };
%-- end accepted_funcs defn.---------
if isstruct(method)
    if length(usave) < usave_length_orig+length(method)
        for i = usave_length_orig+1:usave_length_orig+length(method)
            usave(usave_length_orig+i) = fill_usave(usave(usave_length_orig-1+i));
        end
    end
    for i = 1:length(method)
        if ~strmatch(method(i).func,accepted_funcs)
            error(['Method #' num2str(i) 'is not acceptable.\n  Please use one of the following accepted method.func values:\n' accepted_funcs]);
        end
        usave(usave_length_orig+i).method = method(i);
        if strmatch('write', method(i).func, 'exact')
          usave(usave_length_orig+i).method.nanpix = usave(usave_length_orig).method.nanpix
        end        
    end
else
    error('you must supply a method structure.\n');
    % usave(end) = default(usave(end), 'method', struct('frequency','stim','n_avg',-4,'func','rigid_nonrigid'));
end
[usave.infile] = deal(infile);
[usave.basefilename] = deal(basefilename);
[usave.outfiles] = deal({outfile}); % save the 'usave' outfile in first cell of usave.outfiles
if usave_length_orig > 0
    if isfield(usave(usave_length_orig),'bad_stacks')
        [usave.bad_stacks] = deal(usave(usave_length_orig).bad_stacks);
        [usave.bad_frames] = deal(usave(usave_length_orig).bad_frames);
        [usave.bad_stack_file] = deal(usave(usave_length_orig).bad_stack_file);
%         [usave.is_good] = deal(usave(usave_length_orig).is_good);
%         [usave.fix_range] = deal(usave(usave_length_orig).fix_range);
%         [usave.fix_index] = deal(usave(usave_length_orig).fix_index);       
    end
end

clear temp temp2 basefilename;

%% Interpret options
if ~exist('options','var')
    oplen = 0;
else
    oplen = length(options);
end
metlen = length(usave)-usave_length_orig;
if oplen == metlen
    for i = 1:oplen
        fn = fieldnames(options(i));
        for j = 1:length(fn)
            usave(usave_length_orig+i).(fn{j}) = options(i).(fn{j});
        end
    end
elseif oplen == 0
    % do nothing (use defaults within function call) if oplen = 0;
elseif oplen < metlen
    for i = 1:metlen
        if i <= oplen
            fn = fieldnames(options(i));
            for j = 1:length(fn)
                usave(usave_length_orig+i).(fn{j}) = options(i).(fn{j});
            end
        else
            fn = fieldnames(options(end));
            for j = 1:length(fn)
                usave(usave_length_orig+i).(fn{j}) = options(end).(fn{j});
            end
        end
    end        

else
    error('options struct cannot be longer than method struct\n');
end
%% Call individual method functions
%
% NOTE: methods should be designed with the following syntax:
%       usave = this_method(smm_in, usave_in)
%         where usave is a struct (not struct array)
%         smm_in is the smm object of the stack to be registered
%         usave_in is a pre-organized struct containing the options described above 
%

% convert all usave structures to standard format
clear temp;
% for i = 1:size(usave,2)   
%     temp(i) = fill_usave(usave(i)); % at this point this is redundant, but leaving in as a check on struct integrity
% end
% usave = temp;
% clear temp;
for i = usave_length_orig+1:length(usave)
    thiscall = str2func(usave(i).method.func);
    temp = fill_usave(thiscall(smm, usave(i)));
    usave(i) = temp;
    save(outfile, 'usave', '-mat', '-v7.3');
end
end


%% begin individual function definitions (ADD TO END, please)
function usave_out = fill_usave(usave)
if isstruct(usave)
    temp = usave;
    usave_out = struct;
elseif isempty(usave)
    temp = struct;
    usave_out = struct;
end
usave_out = default(usave_out, 'infile', []);
usave_out = default(usave_out, 'method', []);
usave_out = default(usave_out, 'bad_frames', []);
usave_out = default(usave_out, 'bad_stack_file', []);
usave_out = default(usave_out, 'bad_stacks', []);
usave_out = default(usave_out, 'base_stack', []);
usave_out = default(usave_out, 'lambda', []);
usave_out = default(usave_out, 'multigrid_ops', []);
usave_out = default(usave_out, 'n_cycles', []);
usave_out = default(usave_out, 'orig_err', []);
usave_out = default(usave_out, 'outfiles', []);
usave_out = default(usave_out, 'rmg', []);
usave_out = default(usave_out, 'stacknum', []);
usave_out = default(usave_out, 'u_err', []);
usave_out = default(usave_out, 'u_matrix', []);
usave_out = default(usave_out, 'basefilename', []);
% usave_out = default(usave_out,'is_good', []);
% usave_out = default(usave_out,'fix_range', []);
% usave_out = default(usave_out,'fix_index', []);
fn = fieldnames(usave);
for i = 1:length(fn)
    usave_out.(fn{i}) = usave.(fn{i});
end
% [p, fields, field_shape, sbase] = extract_fields(usave_out, usave);
% sbase = usave;
% usave_out = fill_fields(fields, field_shape, p, sbase);q
end


function usave_out = rigidthenmultigrid(smm_in, usave_in)
% gather stack info
stack_size = smm_in.size;
header_smm = smm_in.header;
stimuli = header_smm.stim_lookup;
Number_of_stacks = stack_size(end);
stimuli = stimuli(1:Number_of_stacks);
% set up padmat
size_smm = double(smm_in.size);
padmat = ones(1, size_smm(2), size_smm(3), 'single');

usave_out = usave_in;
clear usave_in; % done to try to conserve memory

% identify elements of the stacks to be registered
% step 1 extract parameters from 'method' structure
if isfield(usave_out.method,'navg');
    if ~isempty(usave_out.method.navg)
        navg = usave_out.method.navg; % this means the user wants to average a certain # of frames relative to each registration point
    else
        navg = -4;
        usave_out.method.navg = navg;
    end
else
    navg = -4;
    usave_out.method.navg = navg;
end
if isfield(usave_out.method,'frequency')
    if ~isempty(usave_out.method.frequency)
        freq = usave_out.method.frequency;
    else
        freq = 'stim';
        usave_out.method.frequency = freq;
    end
else
    freq = 'stim';
    usave_out.method.frequency = freq;
end
usave_out.method = default(usave_out.method,'skiprigid',false); % if user wants to skip rigid registration step
usave_out.method = default(usave_out.method,'tricknan',false);  % if user wants  to avoid NaNs by decimation via stack padding
usave_out.method = default(usave_out.method,'output','none');   % set to string variable other than 'none' if you want to write an output file
usave_out.method = default(usave_out.method,'max_cycles', 8); % TODO: add max_cycles to acceptable method structure fields
usave_out.method = default(usave_out.method,'exit_percentage',1); % TODO: add exit_percentage to acceptable method structure fields

% if an old u_matrix exists, find the last filled u_matrix and stacknum and copy them
if ~isempty(usave_out.u_matrix)
    old.stacknum = usave_out.stacknum;
    old.u_matrix = usave_out.u_matrix;
    usave_out.stacknum = [];
    usave_out.u_matrix = cell(1);
else
    old = struct;
end

if isfield(old,'u_matrix')
    if isstr(old.u_matrix)
        if ~isempty(strmatch('linear_uinterp_not_saved', old.u_matrix,'exact'))
            % find the u_matrix before the unsaved (interpolated, presumably) u_matrix
            tmp = load(usave_out.outfiles{1});
            for z = 1:length(tmp.usave)-1
                if iscell(tmp.usave(z).u_matrix) && isstr(tmp.usave(z+1).u_matrix)
                    if ~isempty(strmatch('linear_uinterp_not_saved', tmp.usave(z+1).u_matrix,'exact'))
                        old.stacknum = tmp.usave(z).stacknum;
                        old.u_matrix = tmp.usave(z).u_matrix;
                    end
                end
            end
            clear tmp z;
        end
    usave_out.u_matrix = cell(1);
    end
end

% identify the "base stack" to which all others will be registered
if isempty(usave_out.base_stack)
    if ischar(freq)
        if ~isempty(strmatch('stim', freq, 'exact'))
            base_stack = Number_of_stacks/2;
            this = stimuli(base_stack);
            if this ~= 0
                while this ~= 0
                    base_stack = base_stack - 1;
                    this = stimuli(base_stack);
                end
            else
                while this == 0
                    base_stack = base_stack+1;
                    this = stimuli(base_stack);
                end
                base_stack = base_stack - 1;
            end
        elseif ~isempty(strmatch('all',freq,'exact')) || ~isempty(strmatch('full',freq,'exact'))
            base_stack = Number_of_stacks/2;
            this = stimuli(base_stack);
            if this ~= 0
                while this ~= 0
                    base_stack = base_stack - 1;
                    this = stimuli(base_stack);
                end
            else
                while this == 0
                    base_stack = base_stack+1;
                    this = stimuli(base_stack);
                end
                base_stack = base_stack - 1;
            end
        end
    else
        n_divs = floor(Number_of_stacks/freq);
        base_stack = n_divs/2*freq;
        this = stimuli(base_stack);
        while this ~= 0 || stimuli(base_stack+navg) ~= 0
            base_stack = base_stack-1;
            this = stimuli(base_stack);
        end        
    end 
else
    base_stack = usave_out.base_stack;
end
usave_out.base_stack = base_stack;

% identify the stack_registration_index
if ischar(freq)
    if ~isempty(strmatch('stim',freq,'exact'))
        i = base_stack;
        count = 1;
        while i>=1
            if  stimuli(i) ~=0 && stimuli(i-1) == 0
                stack_reg_index(count) = i - 1;
                count = count+1;
            end
            i = i-1;
        end
        i = base_stack;
        while i<=Number_of_stacks
            if  stimuli(i) ~=0 && stimuli(i-1) == 0
                stack_reg_index(count) = i - 1;
                count = count+1;
            end
            i = i+1;
        end
    elseif ~isempty(strmatch('full',freq,'exact')) || ~isempty(strmatch('all',freq,'exact'))
        count = 1;
        for i = base_stack-1:-1:1
            stack_reg_index(count) = i;
            count = count+1;
        end
        stack_reg_index(count) = base_stack;
        count = count+1;
        for i = base_stack+1:Number_of_stacks
            stack_reg_index(count) = i;
            count = count+1;
        end
    else
        error('did not recognize method.frequency argument.\n');
    end
else
    count = 1;
    for i = base_stack-freq:-freq:1
        stack_reg_index(count) = i;
        count = count+1;
    end
    stack_reg_index(count) = base_stack;
    count = count+1;
    for i = base_stack+freq:freq:Number_of_stacks
        stack_reg_index(count) = i;
        count = count+1;
    end
end
usave_out.stacknum = stack_reg_index;

% set up a cell array to hold stacknums to average
for i = 1:length(stack_reg_index)
    if (stack_reg_index(i))+sign(navg)*(abs(navg)-1) < 1
        stack_avg_index{i} = stack_reg_index(i):sign(navg):1;
    elseif (stack_reg_index(i))+sign(navg)*(abs(navg)-1) > Number_of_stacks
        stack_avg_index{i} = stack_reg_index(i):sign(navg):Number_of_stacks;
    else
        stack_avg_index{i} = stack_reg_index(i):sign(navg):(stack_reg_index(i))+sign(navg)*(abs(navg)-1);
    end
end

% set rmg_params
if ~isempty(usave_out.multigrid_ops)
    ops = usave_out.multigrid_ops;
else
    ops = struct('display','none'); % keep display off for default
    fprintf('Note: register_multigrid_options will be called with default values.\n');
    % do not save this fact, just leave it blank%
end
% set lambda
if ~isempty(usave_out.lambda)
    lambda = usave_out.lambda;
else
    lambda = 1e5;
    usave_out.lambda = lambda;
end
% set base image (imBase)
sbi = find(stack_reg_index == base_stack,1);
imBase = double(nanmean(smm_in(:,:,:,stack_avg_index{sbi}),4));

% if trying to trick NaN problems (registration code as of 2010_05_25: JPM)
if isfield(usave_out.method,'tricknan')
    if usave_out.method.tricknan
        tMask = ones(size(imBase));
        if isfield(ops,'options')
            if isfield(ops.options,'zero_nans')
                tmp = ops.options.zero_nans;
            end
        end
        ops.options.zero_nans = false;
        tRmg = register_multigrid_options(imBase,tMask,ops);
        if exist('tmp','var')
            ops.options.zero_nans = tmp;
            clear tmp;
        end
        j = length(tRmg.image_grid);
        while isempty(tRmg.image_grid(j).imFixed)
            j = j-1;
        end
        lastnonnan = j;
        nannedout = zeros(tRmg.image_grid(lastnonnan).sz);
        for i = 1:tRmg.image_grid(lastnonnan).sz(3)
            nannedout(:,:,i) = isnan(tRmg.image_grid(lastnonnan).imFixed(:,:,i));
        end
        % find extent of NaN invasion into image at coarsest u grid
        % (assumes symmetrical frame)
        nanframe(1) = size(find(~nansum(nansum(tRmg.image_grid(lastnonnan).imFixed(:,:,:),3),2)),1)/2;
        nanframe(2) = size(find(~nansum(nansum(tRmg.image_grid(lastnonnan).imFixed(:,:,:),3),1)),2)/2;
        nanframe(3) = size(find(~nansum(nansum(tRmg.image_grid(lastnonnan).imFixed(:,:,:),2),1)),1)/2;
        nanratio(1) = nanframe(1)/tRmg.image_grid(lastnonnan).sz(1);
        nanratio(2) = nanframe(2)/tRmg.image_grid(lastnonnan).sz(2);
        nanratio(3) = nanframe(3)/tRmg.image_grid(lastnonnan).sz(3);
        clear tMask nannedout nanframe;
        % pad the original stack with zeros around the edges to soak up NaN space
        nanpix(1) = ceil(nanratio(1)*tRmg.image_grid(1).sz(1));
        nanpix(2) = ceil(nanratio(2)*tRmg.image_grid(1).sz(2));
        nanpix(3) = ceil(nanratio(3)*tRmg.image_grid(1).sz(3));
        nanpad{1} = zeros([tRmg.image_grid(1).sz(1) nanpix(1) tRmg.image_grid(1).sz(3)]);
        nanpad{2} = zeros([nanpix(2) tRmg.image_grid(1).sz(2)+nanpix(1)*2 tRmg.image_grid(1).sz(3)]);
        nanpad{3} = zeros([tRmg.image_grid(1).sz(1)+nanpix(1)*2 tRmg.image_grid(1).sz(2)+nanpix(2)*2 nanpix(3)]);
        imBase = padnans(imBase, nanpad);
        clear tRmg;
    end
else
    usave_out.method.tricknan = false;
end

% set "real" mask
if ~isfield(usave_out.method,'deadframes')
    if exist('nanpix','var')
        usave_out.method.deadframes = nanpix(3);
    else
        usave_out.method.deadframes = 0;
    end
else 
    if exist('nanpix','var') 
        usave_out.method.deadframes = usave_out.method.deadframes + nanpix(3);
    end
end
if ~isfield(usave_out.method,'deadperim')
    if exist('nanpix','var')
        usave_out.method.deadperim = nanpix(1);
    else
        usave_out.method.deadperim = 0;
    end
else
    if exist('nanpix', 'var')
        usave_out.method.deadperim = usave_out.method.deadperim + nanpix(1);
    end
end
if ~isfield(usave_out.method,'filter')
    usave_out.method.filter = fspecial('gaussian',16,4);
end
% attempt to match masks across improvements
if ~isempty(usave_out.rmg)
    if sum(abs(diff(cat(3,size(usave_out.rmg.image_grid(1).mask),size(imBase)),[],3)))==0
        imMask = usave_out.rmg.image_grid(1).mask;
    else
        warning('A different mask will be applied at this stage compared to previous stages.  You may want to check your stack_automask parameters now. Be careful!\n');
        imMask = stack_automask(imBase,usave_out.method);
    end
else
    imMask = stack_automask(imBase,usave_out.method);
end
usave_out.rmg = register_multigrid_options(imBase,imMask,ops);

% look for evidence of "rnr_temp_u.mat" file (i.e. in-progress)
temp = find(usave_out.outfiles{1}=='/');
outdir = usave_out.outfiles{1}(1:temp(end));
tempfile = [outdir 'rnr_temp_u.mat'];
if exist(tempfile,'file')
    fprintf('\nFound existing "rnr_temp_u.mat" file, continuing previous.  Delete temp file if you wish to start over.\n');
    load(tempfile,'-mat');
    starti = length(usave_out.u_matrix)+1;
else
    starti = 1;
end
clear temp;

% open .cam file for writing if user wants it
writecam = false;
if ~isempty(strmatch('none',usave_out.method.output))
    writecam = false;
elseif ~ischar(usave_out.method.output)
    error('''method.output'' must contain a string variable.\n');
elseif ~isempty(usave_out.basefilename);
  base = usave_out.basefilename;
  camfile = [base '_' usave_out.method.output '.cam'];
  writecam = true;
else
    temp = findstr('.imagine',usave_out.infile);
    if ~isempty(temp)
        base = usave_out.infile(1:temp(end)-1);
        camfile = [base '_' usave_out.method.output '.cam'];
    end
    temp = findstr('.cam', usave_out.infile);
    if ~isempty(temp)
        base = usave_out.infile(1:temp(end)-1);
        camfile = [base '_' usave_out.method.output '.cam'];
    end
    if ~exist('base', 'var')
        base = usave_out.infile;
        camfile = [base '_' usave_out.method.output '.cam'];
    end
    writecam = true;
end
if writecam == true
    if exist(camfile,'file')
        response = input(['File ' camfile ' already exists, overwrite, append, or cancel? (o/a/c): '],'s');
        if ~isempty(strmatch('o',response)) || ~isempty(strmatch('O',response))
            fid = fopen(camfile,'w');
        elseif ~isempty(strmatch('a',response)) || ~isempty(strmatch('A',response))
            fid = fopen(camfile,'a');
        else
            writecam = false;
        end
    else
        fid = fopen(camfile,'w');
    end
end

% cycle!
tic;
if starti == 1
    thisu = [];
    contin = false;
else
    thisu = usave_out.u_matrix{starti-1};
    contin = true;
end
thisdx = [];
fprintf('rigid-nonrigidly registering stacks. Legend:(stack/n_cycles/resid_err)\n..');
finegrid_size = usave_out.rmg.image_grid(usave_out.rmg.gap_data+1).sz;
for i = starti:length(stack_reg_index)
    if i == sbi % no registration for base_stack
        thisu = [];
        thisdx = [];
        usave_out.orig_err(i) = 0;
        usave_out.u_err(i) = 0;
        usave_out.u_matrix{i} = zeros([finegrid_size 3]);
        usave_out.n_cycles(i) = 0;
        tempImg = imBase; %trial!!! remove if buggy  (JPM)
    end
    % make image to be registered
    these_idx = stack_avg_index{i};
    found_bad = intersect(these_idx,usave_out.bad_stacks);
    if isempty(found_bad)
        thisImg = double(nanmean(smm_in(:,:,:,these_idx),4));
        if usave_out.method.tricknan
            thisImg = padnans(thisImg,nanpad);
        end
    else % modify input stacks to incorporate "repaired" bad stacks
      smm_bad = stackmm(usave_out.bad_stack_file);
      [~, bad_in_these, bad_idx] = intersect(these_idx, usave_out.bad_stacks);
      %bad_idx = find(usave_out.bad_stacks == found_bad);
      good_idx = these_idx;
      good_idx(bad_in_these)=[];
      repaired_stacks = double(smm_bad(:,:,:,bad_idx));
      good_stacks = double(smm_in(:,:,:,good_idx));
      thisImg = nanmean(cat(4,repaired_stacks,good_stacks),4);
      if usave_out.method.tricknan
        thisImg = padnans(thisImg,nanpad);
      end
    end
    % if we need to continue in the middle, need to make a copy of the
    % warped lastImg from previous registration try
    if contin == true && ~usave_out.method.skiprigid
        these_idx = stack_avg_index{i-1};
        found_bad = intersect(these_idx,usave_out.bad_stacks);
        if isempty(found_bad)
            tempImg = double(nanmean(smm_in(:,:,:,these_idx),4));
            if usave_out.method.tricknan
                tempImg = padnans(tempImg,nanpad);
            end
        else % modify input stacks to incorporate "repaired" bad stacks
          smm_bad = stackmm(usave_out.bad_stack_file);
          [~, bad_in_these, bad_idx] = intersect(these_idx, usave_out.bad_stacks);
          %bad_idx = find(usave_out.bad_stacks == found_bad);
          good_idx = these_idx;
          good_idx(bad_in_these)=[];
          repaired_stacks = double(smm_bad(:,:,:,bad_idx));
          good_stacks = double(smm_in(:,:,:,good_idx));
          tempImg = nanmean(cat(4,repaired_stacks,good_stacks),4);
          if usave_out.method.tricknan
            tempImg = padnans(tempImg,nanpad);
          end
        end
        tempImg = register_multigrid_warp(tempImg,thisu,usave_out.rmg);
        %clear tempImg;
    end
    % identify original (input) error
    usave_out.orig_err(i) = imcompare(imBase,{thisImg}, imMask);
    
    % if there is an old u_matrix to build from, fill in the u values for this stack
    if isfield(old,'u_matrix')
        if iscell(old.u_matrix)
            [~, closest] = sort(abs(usave_out.stacknum(i)-old.stacknum));
            thesebrack = [old.stacknum(closest(1)) usave_out.stacknum(i) old.stacknum(closest(2))];
            thisu = ucalc(cat(2,old.u_matrix(closest(1)),old.u_matrix(closest(2))),thesebrack);
        end
        clear closest thesebrack;
    end
        
    if ~usave_out.method.skiprigid
        % do the rigid registration (usually very quick);
        if isempty(thisu) && isempty(thisdx)
            [~, thisdx] = register_rigid(imBase,thisImg);
            thisdx = reshape(thisdx,[ones(1,3) 3]);
            thisu = zeros([finegrid_size 3]);
            deci_u = [];
            for j = 1:usave_out.rmg.gap_data+1
                deci_u(j,:) = usave_out.rmg.image_grid(j).restrict;
            end
            deci_u = sum(deci_u,1);
            for j = 1:3
                thisu(:,:,:,j) = thisdx(j)/(2^deci_u(j));
            end
        else
            % register rigid works from previous registered stuff, adds rigidly to previous u values
            if exist('tempImg', 'var');
              [~, thisdx] = register_rigid(tempImg.*imMask,register_multigrid_warp(thisImg,thisu,usave_out.rmg).*imMask);
              %[~, thisdx] = register_rigid(tempImg.*imMask,thisImg.*imMask);              
              thisdx(isnan(thisdx)) = 0;
              thisdx = reshape(thisdx,[ones(1,3) 3]);
              thisu = usave_out.u_matrix{i-1}; % this is the starting spot
              deci_u = [];
              for j = 1:usave_out.rmg.gap_data+1
                deci_u(j,:) = usave_out.rmg.image_grid(j).restrict;
              end
              deci_u = sum(deci_u,1);
              for j = 1:3
                thisu(:,:,:,j) = thisu(:,:,:,j) + thisdx(j)/(2^deci_u(j));
              end
            end
        end
        clear deci_u thisdx;
    end
    
    % do the nonrigid registration
    max_cycles = usave_out.method.max_cycles; % TODO: add max_cycles to acceptable method structure fields
    exit_percentage = usave_out.method.exit_percentage; % TODO: add exit_percentage to acceptable method structure fields
    this = 1;
    cycle = 1;
    err_diff = [];
    pct_err = [];
    while abs(100*this) > exit_percentage && cycle <= max_cycles
        [thisu, err] = register_multigrid_vcycle(thisu,thisImg,lambda,usave_out.rmg);
        tempImg = register_multigrid_warp(thisImg,thisu,usave_out.rmg);
        err = imcompare(imBase, {tempImg}, imMask);
        pct_err(cycle) = err/usave_out.orig_err(i);%err_improve(cycle) = 1-(usave_out.orig_err - err) / usave_out.orig_err;
        if cycle == 1
            err_diff(cycle) = 1-pct_err(cycle);
        else
            err_diff(cycle) = pct_err(cycle-1)-pct_err(cycle);
        end
        if ~isempty(strmatch('all',usave_out.rmg.display,'exact'))
            fprintf('error improvement on full pass #%d is to: %4.1f%%, change of %4.1f%% from previous pass.  ', cycle, pct_err(cycle)*100, 100*err_diff(cycle));
            fprintf('\n');
        end
        this = err_diff(cycle); if isnan(this); this = 1; end; if isinf(this); this = 0; end;
        cycle = cycle+1;
    end
    
    % end extra diagnostic
    usave_out.n_cycles(i) = cycle-1;
    usave_out.u_err(i) = err;
    usave_out.u_matrix{i} = thisu;
    if exist('nanpix','var')
        usave_out.method.nanpix = nanpix;
    end
    
    % write data to .cam file if user says to
    if writecam 
        warpnwrite(fid, smm_in, thisu, usave_out.stacknum(i), usave_out)
    end
    
    % write u-values to a temporary file *'rnr_temp_u.mat'* in same directory as outfile
    save(tempfile, 'usave_out', '-mat', '-v7.3');
    
    % user output
    fprintf('%d/%d/%4.1f..\n', usave_out.stacknum(i),usave_out.n_cycles(i),usave_out.u_err(i)/usave_out.orig_err(i)*100);
    if mod(i,5)
        fprintf('\n');
        toc;
    end
    clear temp;
end
% close any files you may have written to
fclose('all');
% delete tempfile
delete(tempfile);
end


function usave_out = u_interpolate(smm_in, usave_out)
if ~isfield(usave_out.method,'mode')
    usave_out.method.mode = 'linear';
    fprintf('No ''mode'' field input for ''u_interpolate'', using default ''linear''.\n');
end
ops.mode = usave_out.method.mode;
if ~isfield(usave_out.method,'startstack')
    usave_out.method.startstack = 1;
    fprintf('No ''startstack'' field input for ''u_interpolate'', using default (1).\n');
end
ops.startstack = 1;
if ~isfield(usave_out.method,'output')
    usave_out.method.output = 'none';
    if isempty(usave_out.method.output) || ~isempty(strmatch(usave_out.method.output, ' ', 'exact'))
        usave_out.method.output = 'none'; % to make sure a good .cam file isn't overwritten
    end
    fprintf('No ''output'' field input for ''u_interpolate'', using default ''none''.\n');
end

usave_out = uinterp(smm_in, usave_out, ops);

fprintf('Done!\n');
fclose('all');
end

function uout = ucalc(u_in, urange)
% returns the u_value interpolated between two u_in
   if urange(2) == urange(1)
       uout = u_in{1};
   elseif urange(2) == urange(3)
       uout = u_in{2};
   else
       udiff = diff(cat(5,u_in{1},u_in{2}),[],5);
       divisor = abs(urange(end)-urange(1))-1;
       uout = u_in{1}+abs(urange(1)-urange(2))*udiff/divisor;
   end
end

function usave_out = write(smm_in, usave_in)
stack_size = smm_in.size;
header_smm = smm_in.header;
index = usave_in.stacknum;
Number_of_stacks = stack_size(end);

usave_out = usave_in;
clear usave_in; % done to try to conserve memory
% open .cam file for writing if user wants it
writecam = false;
if ~isempty(strmatch('none',usave_out.method.output))
    writecam = false;
elseif ~ischar(usave_out.method.output)
    error('''method.output'' must contain a string variable.\n');
elseif ~isempty(usave_out.basefilename);
  base = usave_out.basefilename;
  camfile = [base '_' usave_out.method.output '.cam'];
  writecam = true;
else
    temp = findstr('.imagine',usave_out.infile);
    if ~isempty(temp)
        base = usave_out.infile(1:temp(end)-1);
        camfile = [base '_' usave_out.method.output '.cam'];
    end
    temp = findstr('.cam', usave_out.infile);
    if ~isempty(temp)
        base = usave_out.infile(1:temp(end)-1);
        camfile = [base '_' usave_out.method.output '.cam'];
    end
    if ~exist('base', 'var')
        base = usave_out.infile;
        camfile = [base '_' usave_out.method.output '.cam'];
    end
    writecam = true;
end
if writecam == true
    if exist(camfile,'file')
        response = input(['File ' camfile ' already exists, overwrite, append, or cancel? (o/a/c): '],'s');
        if ~isempty(strmatch('o',response)) || ~isempty(strmatch('O',response))
            fid = fopen(camfile,'w');
        elseif ~isempty(strmatch('a',response)) || ~isempty(strmatch('A',response))
            fid = fopen(camfile,'a');
        else
            writecam = false;
        end
    else
        fid = fopen(camfile,'w');
    end
    fprintf('Stack: ')
    for i = 1:Number_of_stacks
      fprintf('%d...',i);
      thisu = usave_out.u_matrix{index(i)};
      if writecam 
        warpnwrite(fid, smm_in, thisu, usave_out.stacknum(index(i)), usave_out)      
      else
        error('Insufficient parameters to write cam file.\n');
      end
    end    
end

 
end
