% Registering stacks

% 2006-09-30

% This script will register stacks using multires with given g0
% multires with rigid transformation for stimuli data

% copywrite Diwakar Turaga 2007

clear all

%--------------------------------------------------------------------------
% INPUT PARAMETERS ... make sure they are all looking right

datadir = '/mnt/dturaga_004/dturaga/2007_10_08/';
savedir = '/mnt/dturaga_004/dturaga/2007_10_08/';
datafile = '2007_10_08_1';
savefile = [datafile '_rigid_reg.cam'];

base_stack_index = 255; % can specify base_stack here, if '0' then will 
                      % automatically calculate base stack somewhere in 
                      % the middle of the stacks
start_stack_number = 1; % the first stack is usually bad looking. so we 
                        % start from stack #2
stop_stack_number = 530;                   

%------------------------------------------------------------------------
% Making sure the files are not over-written

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

size_smm = double(smm.size);
padmat = ones(1, size_smm(2), size_smm(3), 'single');
% padding because stackmm removes one row


%--------------------------------------------------------------------------
% computing the base stack
stack_base = smm(:,:,:,base_stack_index);
stack_base = single(stack_base);
stack_base = imfilter_gaussian(stack_base, [6 6 1]);

%--------------------------------------------------------------------------
% Warping the stacks

fprintf('registering stacks\n');


for i = start_stack_number:stop_stack_number
    
    stack_temp= single(smm(:,:,:,i));
    stack_temp = imfilter_gaussian(stack_temp, [6 6 1]);
    
    options.pixel_spacing = [0.71 0.71 6];
    options.n_vcycles = 4;
    %options.display = 1;
    
    [params,im] = register_translation_multires(stack_base,stack_temp, options);
    %[img, params] = register_rigid(stack_base, stack_temp);
    
    params
    
    img = single(smm(:,:,:,i));
    img = image_shift(img, params);

    
    img_cat = cat(1, img, padmat);
    count = fwrite(fid,uint16(img_cat),'uint16');
    clear stack_temp_warped img img_cat
    
    fprintf('..%d..',i);
end
fprintf('\n')
fprintf('DONE !\n')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------