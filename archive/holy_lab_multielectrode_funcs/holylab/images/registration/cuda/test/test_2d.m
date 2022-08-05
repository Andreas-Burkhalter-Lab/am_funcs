clear;
cuda_cleanup;
close all;
clc;

%% Create the fixed & moving 2D images
% im = double(imread('cameraman.tif'));
% scale = 6;  % max # of pixels in shift
% sz = size(im);
% udef = randn([sz 2]);
% udef = imfilter_gaussian(udef,[60 60 0]);  % smooth the random deformation
% udef = udef * (scale / max(abs(udef(:)))); % set the scale
% imw_full = register_multigrid_warp(im,udef);
% rng = 15:240;
% moving = imw_full(rng,rng);
% fixed = im(rng,rng);

%% or load the saved fixed & moving 2D images
load /home/jian/Documents/sample_data/data61.mat;

%% RGB compare the fixed & moving 2D images
figure; imshowrgb(fixed,moving); title('Unregistered fixed & moving images');

%% Set the registration parameters
pyramid = array_restrict_schedule(size(fixed));
ops = register_phasecorr_initialize(fixed,struct('pyramid',pyramid));
ops.lambda = 1.0E+1; % suggested value
ops.plot = false;
ops.barrier = true;
ops.barrier_mu = 1.0E-9; % suggested value

%% Repeat the block improve registration using a coarse to fine grid
length_pyramid = prod(size(pyramid));
for i = length_pyramid:-1:(length_pyramid-4)
  % Register the image
  usz = pyramid(i).sz;
  if int8(i) == int8(length_pyramid)
    uc = register_block_improve(usz,fixed,moving,ops);
  else
    uc{1} = array_prolong(uc{1}, usz);
    uc{2} = array_prolong(uc{2}, usz);
    uc = register_block_improve(uc,fixed,moving,ops);
  end
  % RGB compare the fixed & moving 2D images
  imw = register_phasecorr_warp(uc,moving,ops);
  j = 0;
  figure; imshowrgb(fixed,imw); title(sprintf('[%d %d %d]',size(uc{1}), j));
  
  continue;
end

% Also show the deformation itself
% funch = register_gui_utilities;
% unorm = funch.upc2mg(uc,ops.sz_spatial);
% figure; imflow_display(unorm);

