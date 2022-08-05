% This is a short demo script on running multigrid registration.  There
% is no help for this, you simply need to look at the code.

% Get the fixed and moving images
smm = stackmm('vno_2009_12_29_2.imagine');
imM = double(smm(:,:,:,500));  % moving image
imF = double(smm(:,:,:,1));  % fixed image
% We are going to ignore pixels that lie "near" the edges, so define what
% "near" means along each coordinate
trimsz = [50 50 4];
mask = imF > 1500;  % 1500 here represents a cutoff for "black" pixels, these
                    % are useless so don't waste computing time on these.
mask = trim_mask_edges(mask,trimsz);
% lambda controls how much you are willing to tolerate volume changes; it
% is the coefficient of the volume-change penalty (see
% register_logdetpenalty).  A good initial guess would be something on
% the order of mean(imF(:)^2).  Smaller values make you more tolerant of
% volume changes, larger values force a more rigid-like deformation.
lambda = 1e5;
% Now make the grid of coarser images
rmg_params = register_multigrid_options(imF,mask,struct('pixel_spacing',[0.71 0.71 200/39]));
rmg_params.gap_data = 4;  % this is an important parameter, see the help
% The next is optional, but it can help...Start with a _rigid_
% registration to get the offset (for now, this can only handle
% translations, not rotations)
[imr,dx] = register_rigid(imF,imM);
u0 = register_rigid2nonrigid(dx,rmg_params);
% If you don't do rigid registration first, just set u0 = [] for what follows
% Now the big calculation!
[u,err,rmgp] = register_multigrid_vcycle(u0,imM,lambda,rmg_params);
% Evaluate the warped image, so you can see it.
imMg = register_multigrid_warp(imM,u,rmg_params);
% Trim off the edges, since we didn't try to match them anyway.
imFmask = imF .* mask;
imMgmask = imMg .* mask;
% Look at the result
figure; for i = trimsz(3)+1:40-trimsz(3)
  imshowrgb(imFmask(:,:,i),imMgmask(:,:,i));
  title([num2str(i) '/40']);
  pause
end
% You can look at the deformation itself by calling register_visualize_u.
