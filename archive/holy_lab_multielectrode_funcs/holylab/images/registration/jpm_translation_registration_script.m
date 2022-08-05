% Registering stacks
% 2010-01-27

% This script will register stacks using translation with given dx

% Adapted from original DT script by JPM


%--------------------------------------------------------------------------
%% INPUT PARAMETERS ... make sure they are all looking right

datadir = '/mnt/hammeng_01b1/hammeng/OCPI/aob_img/2010_01_25/';
savedir = '/home/julian/';
datafile = '2010_01_25_GAD65_prep1_anterior_timetest';
savefile = [datafile '_reg.cam'];
savedir_dx = [savedir 'dx'];
savedir_dx_linear = [savedir 'dx_linear_' datafile];


base_stack_index = 450; % can specify base_stack here, if '0' then will 
                      % automatically calculate base stack somewhere in 
                      % the middle of the stacks
start_stack_number = 2; % the first stack is usually bad looking. so we 
                        % start from stack #2
bad_stacks = [239 541 674];  % ignore these stacks (dropped frames, etc)
stop_stack_number = 900;
reg_wrt_stimuli = 1; % this is the number of stack before the stimuli 
                     % you want to register, normally one
stack_spacing = 10;  % number of stacks between registration attempts
options.pixel_spacing = [1 1 10];
mask = 0; % This will create a mask

% ALSO LOOK AT THE MASK keep_range = [20 -50]; This looks at height above
% and below the surface of the tissue the registration is done.

%------------------------------------------------------------------------
%% Making sure the .cam file is not over-written
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
%% this is to fwrite the stacks before start_stack_number, i.e. usually the
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
%% calculating base_stack_index
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
%% computing the base stack
stack_base = smm(:,:,:,base_stack_index);
stack_base = single(stack_base);

% -------------------------------------------------------------------------
%% Creating mask
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
%% calculating the indeces of stacks which are going to be registered
c = 1;
i = base_stack_index;

if ~isempty(stimuli)
    while i>=start_stack_number
       if  stimuli(i) ~=0 && stimuli(i-1) == 0
           stack_reg_index(c) = i - reg_wrt_stimuli;
           c = c+1;
       end
       i = i-1;
    end
else
    stack_reg_index = base_stack_index-1:-1*stack_spacing:start_stack_number;
    if stack_reg_index(end)~=start_stack_number
        stack_reg_index(end+1) = start_stack_number;
    end
    e1 = size(stack_reg_index,2);
    e2 = size([base_stack_index+1:10:Number_of_stacks],2);
    stack_reg_index(e1+1:e1+e2) = base_stack_index+1:10:Number_of_stacks;
    if stack_reg_index(end)~=Number_of_stacks
        stack_reg_index(end+1) = Number_of_stacks;
    end
    rmstack = indexainb(bad_stacks,stack_reg_index);
    stack_reg_index(rmstack) = [];
end
% 
% stack_reg_index(c) = start_stack_number;
% c = c+1;
% 
% i = base_stack_index;
% while i<=stop_stack_number
%     
%    if  stimuli(i) ~=0 && stimuli(i-1) == 0
%         stack_reg_index(c) = i - reg_wrt_stimuli;
%         c = c+1;
%    end
%    i = i+1;
% end
% 
% stack_reg_index(c) = stop_stack_number;
stack_reg_index % let user see the stacks which are going to be registered

%---------------------------------------------------------------------
%% registering stacks
fprintf('Progress (out of %d stacks):\n',stop_stack_number);

options.n_vcycles = 2;
options.min_pixels = 10;
options.display = false;

for stack_index= 1:length(stack_reg_index)

    stack_temp = single(smm(:,:,:,stack_reg_index(stack_index)));

    if stack_index == 1 || abs(stack_reg_index(stack_index)-base_stack_index)<=1
        options.dx_start = [0 0 0];%register_g0(size(g{1}));
    else
        options.dx_start = dx;
    end
    
    [dx, reg_stack] = register_translation_multires(stack_base, stack_temp, options); 

    save_file_dx = ['dx_' num2str(stack_reg_index(stack_index))];
    if ~exist(savedir_dx, 'dir')
        mkdir(savedir_dx);
    end
    save([savedir_dx filesep save_file_dx], 'dx')

    fprintf('..%d..',stack_reg_index(stack_index));
end
%--------------------------------------------------------------------------
%% Load dx vals & if you feel like it, look at the chosen dx values
display_dx = true; % change to false to suppress plotting

for i = 1:size(stack_reg_index,2)
    temp = load([savedir_dx filesep 'dx_' num2str(stack_reg_index(i)) '.mat']);
    dxvals(stack_reg_index(i),:) = temp.dx;
end

if display_dx == true
    figure(101); delete(get(gcf,'children'));
    plot3(dxvals(stack_reg_index,1), dxvals(stack_reg_index,2), dxvals(stack_reg_index,3),'ko');
    set(get(gca,'xlabel'),'string','x-pixels');
    set(get(gca,'ylabel'),'string','y-pixels');
    set(get(gca,'zlabel'),'string','z-pixels');
end

%------------------------------------------------------------------------
%% linearizing dx values between the calculated stacks

fprintf('Linearizing:  ')

% interpolate dx vals between calculated stacks
for i = start_stack_number:size(dxvals,1)-1
    tempi = 1;
    temps = 0;
    while ~temps
        temps = sum(dxvals(i+tempi,:) == [0 0 0]) ~= 3;
        tempi = tempi+1;
    end
    nexti = i+tempi-1;
    if sum(dxvals(i,:) == [0 0 0]) == 3
        dxvals(i,:) = dxvals(lasti,:)+...
                      (dxvals(nexti,:)-dxvals(lasti,:))*(i-lasti)/(nexti-lasti);
    else
        lasti = i;
    end
end 

fprintf('done\n');

%% the middle by base_stack_index, i sort them out to get them in order for linearizing
stimuli_index = sort(stack_reg_index); % since all stimuli index is broken up in ...

Number_of_stimuli = length(stimuli_index);
for i = 1:Number_of_stimuli-1   

    stimuli_start = stimuli_index(i);
    stimuli_end = stimuli_index(i+1);
    
    dx_start_file = ['dx_' num2str(stimuli_start)];
    dx_end_file = ['dxg_' num2str(stimuli_end)];

    load([savedir_dx filesep dx_start_file]);
    dx_start = dx;
    
    load([savedir_dx filesep dx_end_file]);
    dx_end = dx;
    dx_new = dx;
    dx_cell_size = size(dx_new);
    
                
    for dx_counter = stimuli_start:stimuli_end
        alpha = ((dx_counter - stimuli_start)/(stimuli_end - stimuli_start));
            
        for dx_cell_index = 1:dx_cell_size(2)
            dx_new{dx_cell_index} = ((1 - alpha).*dx_start{dx_cell_index}) + (alpha.*dx_end{dx_cell_index});
        end
            
        save_file_dx = ['dx_' num2str(dx_counter)];
        save([savedir_dx_linear filesep save_file_dx], 'dx_new')
             
    end
end

%--------------------------------------------------------------------------
%% Moving the stacks

fprintf('shifting stacks\n');

for i = start_stack_number:stop_stack_number
    stack_temp= single(smm(:,:,:,i));
    
    img = image_shift(stack_temp, dxvals(i,:));
    
    img_cat = cat(1, img, padmat);
    count = fwrite(fid,uint16(img_cat),'uint16');
    if (count < prod(stack_size(1:end-1)))
        error('Wrote fewer pixels than expected, quitting');
    end
    clear img img_cat
    fprintf('..%d..',i);
end
fprintf('\n')
fprintf('DONE !\n')%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------