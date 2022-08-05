clear;
cuda_cleanup;
close all;
clc;

%% profiling on
profile off;
profile clear;
profile on;

%% Load the image stack and prepare the fixed & moving images
% load /home/jian/Documents/sample_data/20120217_WFUrine_male_la_forregistration.mat
load /home/jian/Documents/sample_data/data1.mat
fixed_res = double(fixed);
moving_res = double(moving);

%% RGB compare each frame of two image stacks before registration
% i = 1;
% while ~isnan(i)
%   imshowrgb(uint16(fixed_res(:,:,i)), uint16(moving_res(:,:,i)));
%   title(num2str(i));
%   i = keystepper(1:40,i);
% end

%% Prepare & set the parameters
pyramid = array_restrict_schedule(size(fixed_res));
ops = register_phasecorr_initialize(fixed_res,struct('pyramid',pyramid));
ops.lambda = 1.0E+3; % Best default value for det(J) penalty
ops.plot = false;
ops.barrier = true;
ops.barrier_mu = 1.0E-4; % Best default value for Edge Barrier penalty

%% Do the image registration using a coarse to fine grid
length_pyramid = prod(size(pyramid));
for i = length_pyramid:-1:(length_pyramid-2)
  % Register the image
  usz = pyramid(i).sz;
  if int8(i) == int8(length_pyramid)
    uc = register_block_improve(usz,fixed_res,moving_res,ops);
  else
    uc{1} = array_prolong(uc{1}, usz);
    uc{2} = array_prolong(uc{2}, usz);
    uc{3} = array_prolong(uc{3}, usz);
    uc = register_block_improve(uc,fixed_res,moving_res,ops);
  end
  
  % RGB compare each frame of two image stacks before registration
%   imw = register_phasecorr_warp(uc,moving_res,ops);
%   i = 1;
%   while ~isnan(i)
%     imshowrgb(uint16(fixed_res(:,:,i)), uint16(imw(:,:,i)));
%     title(num2str(i));
%     i = keystepper(1:40,i);
%   end 
end

%% Save profile
profile off;
profsave;

