% Registering stacks

% 2006-09-30

% This script will register stacks using multires with given g0
% multires with rigid transformation for stimuli data


clear all

%--------------------------------------------------------------------------
% INPUT PARAMETERS ... make sure they are all looking right

datadir = '/sata3ntfs/2006-03-27/';
savedir = '/sata2ext3/dturaga/2006-03-27/';
datafile = '2006-03-27-1';
savefile = [datafile '_reg.cam'];
savedir_g = [savedir 'g_' datafile];
savedir_g_linear = [savedir 'g_linear_' datafile];

options.covariant = true;
lambda = 0; % look at register_multires for its description
base_stack_index = 425; % can specify base_stack here, if '0' then will 
                      % automatically calculate base stack somewhere in 
                      % the middle of the stacks
start_stack_number = 2; % the first stack is usually bad looking. so we 
                        % start from stack #2
stop_stack_number = 900;
reg_wrt_stimuli = 1; % this is the number of stack before the stimuli 
                     % you want to register, normally one                     
base_decimation_schedule = [2 2 1]; % Can't register full size 
                                    % 1000by1000by40 -> 500by500by40
options.pyramid_schedule = [2 2 1; 2 2 1;  3 3 1; 2 2 2];
options.n_levels = 4; % ideally the final calculated g should be 
                      %~20by20by20. so after imreduce by 
                      % base_decimation_schedule and
                      % options.pyramid_schedule, the image should be
                      % 20by20by20. n_levels should length of
                      % pyramid_schedule
mask = 1; % This will create a mask

% ALSO LOOK AT THE MASK keep_range = [20 -50]; This looks at height above
% and below the surface of the tissue the registration is done.

%------------------------------------------------------------------------
% Making sure the files are not over-written
if exist([savedir_g],'file')
    response = input('g_value file exists. Replace? (r/n) [n]:', 's');
    if strcmp(lower(response),'r')
        rmdir([savedir_g],'s');
        rmdir([savedir_g_linear],'s');
        mkdir ([savedir_g]);
        mkdir ([savedir_g_linear]);
    end
else
    mkdir ([savedir_g]);
    mkdir ([savedir_g_linear]);
end

if exist([savedir savefile],'file')
    response = input('.cam File exists. Replace or append? (r/a) [a]:', 's');
    if strcmp(lower(response),'r')
        delete([savedir savefile]);
    end
end

[fid,msg] = fopen([savedir savefile],'a');
if (fid < 0)
    fprintf('%s\n',msg);
    error(['Can''t open file ' savedir savefile ' for writing. Do you have permission?']);
end
%------------------------------------------------------------------------
  
smm = stackmm([datadir datafile]);
stack_size = smm.size;
header_smm = smm.header;
stimuli = header_smm.stim_lookup;
Number_of_stacks = stack_size(end);

size_smm = double(smm.size);
padmat = ones(1, size_smm(2), size_smm(3), 'single');
% padding because stackmm removes one row

%------------------------------------------------------------------------
% this is to fwrite the stacks before start_stack_number, i.e. usually the
% first stack which is usually bad
for base_index = 1:start_stack_number-1
    stack_before_start = single(smm(:,:,:,base_index));
    stack_before_start_cat = cat(1, stack_before_start, padmat);
    count = fwrite(fid,uint16(stack_before_start_cat),'uint16');
    if (count < prod(stack_size(1:end-1)))
        error('Wrote fewer pixels than expected, quitting');
    end
end
clear stack_before_start stack_before_start_cat

%--------------------------------------------------------------------------
% calculating base_stack_index
if base_stack_index == 0
    base_stack_index_gross = round(stop_stack_number/2);
    while base_stack_index_gross >= start_stack_number
    
        if stimuli(base_stack_index_gross) ~= 0 && ...
           stimuli(base_stack_index_gross-1) == 0
            
            base_stack_index = base_stack_index_gross - reg_wrt_stimuli;
            break;
        end
    base_stack_index_gross = base_stack_index_gross - 1;   
    end
end
%--------------------------------------------------------------------------
% computing the base stack
stack_base = smm(:,:,:,base_stack_index);
stack_base = single(stack_base);

for decIndex = 1:size(base_decimation_schedule,1)
    stack_base = imreduce(stack_base, base_decimation_schedule(decIndex,:));
end

% -------------------------------------------------------------------------
% Creating mask
if mask              % if mask = 0, then stack_base is the entire stack_base
    dimorder = [2 1 3];
    stack_base_p = permute(stack_base,dimorder);
    z = tissueboundary_initialize(stack_base_p,1200);
    zmf = medfilt2(z,[11 3]);
    stack_base_mask = false(size(stack_base));
    keep_range = [20 -50]; % 20 pixels above & 30 pixels below
    kr = min(keep_range):max(keep_range);
    for i = 1:size(zmf,1)
        for j = 1:size(zmf,2)
            if (zmf(i,j)+kr(1) > 1 && zmf(i,j)+kr(end) < size(stack_base_p,1))
                stack_base_mask(i,zmf(i,j) + kr,j) = true;
            end
        end
    end
    stack_base_final_mask = stack_base;
    stack_base_final_mask(~stack_base_mask) = nan;
    stack_base = stack_base_final_mask; 
    clear stack_base_final_mask stack_base_mask 
    clear stack_base_p z zmf keep_range kr
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% calculating the indeces of stacks which are going to be registered
c = 1;
i = base_stack_index;

while i>=start_stack_number
   if  stimuli(i) ~=0 && stimuli(i-1) == 0
        stack_reg_index(c) = i - reg_wrt_stimuli;
        c = c+1;
   end
   i = i-1;
end

stack_reg_index(c) = start_stack_number;
c = c+1;

i = base_stack_index;
while i<=stop_stack_number
    
   if  stimuli(i) ~=0 && stimuli(i-1) == 0
        stack_reg_index(c) = i - reg_wrt_stimuli;
        c = c+1;
   end
   i = i+1;
end

stack_reg_index(c) = stop_stack_number;
stack_reg_index % let user see the stacks which are going to be registered

%---------------------------------------------------------------------
% registering stacks
fprintf('Progress (out of %d stacks):\n',stop_stack_number);

for stack_index= 1:length(stack_reg_index)

    stack_temp = single(smm(:,:,:,stack_reg_index(stack_index)));

    for decIndex = 1:size(base_decimation_schedule,1)
        stack_temp = imreduce(stack_temp, base_decimation_schedule(decIndex,:));
    end

    if stack_reg_index(stack_index) == base_stack_index
        options.g0 = register_g0(size(g{1}));
    end
    
    g = register_multires(stack_base, stack_temp,lambda,options);

    options.g0 = g;

    save_file_g = ['g_' num2str(stack_reg_index(stack_index))];
    save([savedir_g filesep save_file_g], 'g')

    fprintf('..%d..',stack_reg_index(stack_index));
end

%------------------------------------------------------------------------
% starting linearizing g values between the calculated stacks

fprintf('Linearizing:\n')

stimuli_index = sort(stack_reg_index); % since all stimuli index is broken up in ...
% the middle by base_stack_index, i sort them out to get them in order for linearizing

Number_of_stimuli = length(stimuli_index);
for i = 1:Number_of_stimuli-1   

    stimuli_start = stimuli_index(i);
    stimuli_end = stimuli_index(i+1);
    
    g_start_file = ['g_' num2str(stimuli_start)];
    g_end_file = ['g_' num2str(stimuli_end)];

    load([savedir_g filesep g_start_file]);
    g_start = g;
    
    load([savedir_g filesep g_end_file]);
    g_end = g;
    g_new = g;
    g_cell_size = size(g_new);
    
                
    for g_counter = stimuli_start:stimuli_end
        alpha = ((g_counter - stimuli_start)/(stimuli_end - stimuli_start));
            
        for g_cell_index = 1:g_cell_size(2)
            g_new{g_cell_index} = ((1 - alpha).*g_start{g_cell_index}) + (alpha.*g_end{g_cell_index});
        end
            
        save_file_g = ['g_' num2str(g_counter)];
        save([savedir_g_linear filesep save_file_g], 'g_new')
             
    end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Warping the stacks

fprintf('warping stacks\n');

stack_base = smm(:,:,:,base_stack_index);
stack_base = single(stack_base);

for i = start_stack_number:stop_stack_number
    
    stack_temp= single(smm(:,:,:,i));
    
    g_index = ['g_' num2str(i)];
    load([savedir_g_linear filesep g_index '.mat'])

    g_hires = register_expandg(g_new,size(stack_temp));
    stack_temp_warped = register_warp(stack_temp,g_hires); 
    clear stack_temp g_hires

    % this part will rigid register the warped stack whenever the stimulus
    % is on. this is to take care of the mechanical jiggering occuring when
    % a valve is opened
    
    
    if stimuli(i) ~= 0
        [img, params] = register_rigid(stack_base, stack_temp_warped);
    else
        img = stack_temp_warped;
    end
    
    img_cat = cat(1, img, padmat);
    count = fwrite(fid,uint16(img_cat),'uint16');
    if (count < prod(stack_size(1:end-1)))
        error('Wrote fewer pixels than expected, quitting');
    end
    clear stack_temp_warped img img_cat
    
    fprintf('..%d..',i);
end
fprintf('\n')
fprintf('DONE !\n')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------