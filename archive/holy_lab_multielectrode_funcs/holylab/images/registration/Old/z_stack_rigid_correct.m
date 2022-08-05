function registered_img = z_stack_rigid_correct(img, varargin)
% function registered_img = z_stack_rigid_correct(img) uses register_rigid
% to neutralize x-y shifts during acquisition
%
% this function takes z-stack images of 3-4 dimensions:
%            img(x_pixels, y_pixels, z_frames, n_channels (optional))
% 
% See also REGISTER_RIGID

% Copyright 2007 Julian P. Meeks

% Interpret and screen img_dimensions (either 3 or 4)

if size(varargin,2) > 0
    if varargin{1} == 'stackmm'
        stackmm = 1;
        img_size = img.size;
        img_dimensions = size(img_size,2);
    end
else
    stackmm = 0;
    img_size = size(img);
    img_dimensions = size(img_size,2);
end

if ~ismember(img_dimensions,[3 4]);
    error('Images must be 3 or 4-dimensional');
    return;
end

% set up order for aligning (starting with center frame)
if stackmm
    n_frames = img_size(3);
else
    n_frames = img_size(3);
end
center_frame = ceil(n_frames/2);

% cycle through frames and using register_rigid
toggle = 1;
pos_indices = center_frame:n_frames;
neg_indices = rot90(1:center_frame,2);
if stackmm &&  img_dimensions == 3
    registered_img = single(img(:,:,:));
elseif stackmm && img_dimensions == 4
    registered_img = single(img(:,:,:,:));
else
    registered_img = single(zeros(size(img)));
end

neg.dx_start = single([0 0]);
pos.dx_start = single([0 0]);

%options=struct('subpixel',true);

if img_dimensions == 3
    registered_img(:,:,center_frame) = single(img(:,:,center_frame));
    for idx_cycle = 2:max([length(pos_indices), length(neg_indices)])
        if idx_cycle <= length(neg_indices)
            [neg.dx_start registered_img(:,:,neg_indices(idx_cycle))] = ...
                register_translation_multires(single(registered_img(:,:,center_frame)),...
                               single(img(:,:,neg_indices(idx_cycle))),neg.dx_start);
        end
        if idx_cycle <= length(pos_indices)
            [pos.dx_start registered_img(:,:,pos_indices(idx_cycle))] = ...
                register_translation_multires(single(registered_img(:,:,center_frame)), ...
                               single(img(:,:,pos_indices(idx_cycle))),pos.dx_start);
        end
        fprintf('..%d',idx_cycle); if mod(idx_cycle, 10)==0; fprintf('\n');end;
    end
elseif img_dimensions == 4
    registered_img(:,:,center_frame,:) = img(:,:,center_frame,:);
    for idx_cycle = 2:max([length(pos_indices), length(neg_indices)])
        if idx_cycle <= length(neg_indices)
            [neg.dx_start registered_img(:,:,neg_indices(idx_cycle),:)] = ...
                register_translation_multires(single(registered_img(:,:,center_frame,:)),...
                               single(img(:,:,neg_indices(idx_cycle),:)),neg.dx_start);
        end
        if idx_cycle <= length(pos_indices)
            [pos.dx_start registered_img(:,:,pos_indices(idx_cycle),:)] = ...
                register_translation_multires(single(registered_img(:,:,center_frame,:)),... 
                               single(img(:,:,pos_indices(idx_cycle),:)),pos.dx_start);
        end
        fprintf('..%d',idx_cycle); if mod(idx_cycle, 10)==0; fprintf('\n');end;
    end
end

end