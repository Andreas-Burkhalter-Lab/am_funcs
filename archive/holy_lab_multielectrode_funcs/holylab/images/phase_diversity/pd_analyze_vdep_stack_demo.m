% pth = '/mnt/dturaga_005/dturaga/2009_05_29/2009_05_29_beads_dm_vol_series_5act/
actuator = 46;
x0 = {605, 718};
snipsize = 63;

% Set up the pupil
NA = 0.5;
M = 200/9;
lambda = 0.5;
pixel_spacing = 6.45;

snipindx = {(-snipsize:snipsize+1)+x0{1},(-snipsize:snipsize+1)+x0{2}};
imsz = 2*(snipsize+1)*[1 1];

% Load the data and cut out the snippet
if (~exist('imsnip','var') || actuator ~= actuatorOld)
  actuatorOld = actuator;
  [pupildata.H0,pupildata.rho,pupildata.theta] = pupil_initialize(NA,M,lambda,imsz,pixel_spacing);
  % Load the aberrated data
  imindex = 1:21;
  Khalf = length(imindex);
  load(['cam1_act' num2str(actuator)])
  imsnip = double(im1stack(snipindx{:},imindex));
  clear im1stack
  load vol_tune
  vol_tune = vol_tune(imindex);
  % also load the unaberrated image (used for registration)
  load(['cam2_act' num2str(actuator)])
  imsnip_unab = double(im2stack(snipindx{:},imindex));
  clear im2stack
  params.im_unab = imsnip_unab;
end

% Show the user the stack
hfig = figure; 
k = 1;
while ~isnan(k)
  imshowsc(imsnip(:,:,k));
  title(sprintf('Frame %d (q to continue)',k))
  k = keystepper(1:size(imsnip,3),k);
end
close(hfig); drawnow

[p,imc] = pd_analyze_vdep_stack(imsnip,vol_tune,pupildata,params);

figure
k = 1;
while ~isnan(k)
  subplot(1,2,1)
  imshowsc(imsnip(:,:,p.selIndex(k)))
  title(sprintf('Original (index %d)',k))
  subplot(1,2,2)
  imshowsc(imc(:,:,k))
  title('Calculated')
  k = keystepper(1:length(p.selIndex),k);
end

