clear;
close all;
clc;

%% Load the image stack and prepare the fixed & moving images
% load /home/jian/Documents/sample_data/20120217_WFUrine_male_la_forregistration.mat
% fixed = stk{30};
% moving = stk{3};
load /home/jian/Documents/sample_data/data1.mat
fixed_res = double(fixed);
moving_res = double(moving);

%% RGB compare each frame of two image stacks before registration
i = 1;
while ~isnan(i)
  imshowrgb(uint16(fixed_res(:,:,i)), uint16(moving_res(:,:,i)));
  title(num2str(i));
  i = keystepper(1:40,i);
end

%% Prepare & set the parameters
pyramid = array_restrict_schedule(size(fixed_res));
ops = register_phasecorr_initialize(fixed_res,struct('pyramid',pyramid));
ops.lambda = 1.0E-6; % Best default value for det(J) penalty
ratio = 10; % The ratio of decreasing ops.lambda parameter
ops.plot = false;
ops.barrier = true;
ops.barrier_mu = 1.0E-4; % Best default value for Edge Barrier penalty

%% Do the image registration using a coarse to fine grid
length_pyramid = prod(size(pyramid));
for i = length_pyramid:-1:(length_pyramid-3)
  % Register the image
  usz = pyramid(i).sz;
  ops.lambda = ops.lambda / ratio;
  if int8(i) == int8(length_pyramid)
    uc = register_block_improve(usz,fixed_res,moving_res,ops);
  else
    uc{1} = array_prolong(uc{1}, usz);
    uc{2} = array_prolong(uc{2}, usz);
    uc{3} = array_prolong(uc{3}, usz);
    uc = register_block_improve(uc,fixed_res,moving_res,ops);
  end
  
  % RGB compare each frame of two image stacks before registration
  imw = register_phasecorr_warp(uc,moving_res,ops);
  i = 1;
  while ~isnan(i)
    imshowrgb(uint16(fixed_res(:,:,i)), uint16(imw(:,:,i)));
    title(num2str(i));
    i = keystepper(1:40,i);
  end
  
  % Also show the deformation itself
% funch = register_gui_utilities;
% unorm = funch.upc2mg(uc,ops.sz_spatial);
% figure; imflow_display(unorm);
  continue;
  
end

