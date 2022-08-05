clear;
close all;

%% load image stack, prepare fixed & moving images
%load /home/tim/Documents/Sample/20120217_WFUrine_male_la_forregistration.mat
load /home/jian/Documents/sample_data/20120217_WFUrine_male_la_forregistration.mat
fixed_res = double(stk{1,30});
stk_size = length(stk);
stk_reg = cell(1,stk_size);

for stk_index = 1:stk_size
    moving_res = double(stk{1,stk_index});
    
    %% rgb compare each frame of two image stacks
    i = NaN;
    while ~isnan(i)
        imshowrgb(uint16(fixed_res(:,:,i)), uint16(moving_res(:,:,i)));
        title(num2str(i));
        i=keystepper(1:40,i);
    end 
    
    %% apply block mismatch
    pyramid = array_restrict_schedule(size(fixed_res));
    ops = register_phasecorr_initialize(fixed_res,struct('pyramid',pyramid));
    usz = [17, 17, 3]; % pyramid restriction size
    
    for reg_index = 1:2
        mismatch = register_block_mismatch_cuda(usz,fixed_res,moving_res,ops); % block mismatch registration
        
        %% compute the minimum in each block
        uc = cell(1,ndims(fixed_res));
        for i = 1:length(uc)
            uc{i} = zeros(usz);
        end
        for blockIndex = 1:numel(mismatch)
            this_sz = size(mismatch{blockIndex});
            this_center = ceil(this_sz/2);
            [~,minIndex] = min(mismatch{blockIndex}(:));
            coords = ind2sub_matrix(this_sz,minIndex);
            for dimIndex = 1:length(this_sz)
                uc{dimIndex}(blockIndex) = coords(dimIndex) - this_center(dimIndex);
            end
        end
        
        %% warp moving image
        imw = register_phasecorr_warp(uc,moving_res,ops);
        
        %% rgb compare each frame of two image stacks
        i = NaN;
        while ~isnan(i)
            imshowrgb(uint16(fixed_res(:,:,i)), uint16(imw(:,:,i)));
            title(num2str(i));
            i=keystepper(1:40,i);
        end
        
        moving_res = imw;
    end
    
    stk_reg{1,stk_index} = single(moving_res);
end

%% save registered image stack
save stk_reg.mat h stacklist stk_reg -v7.3

