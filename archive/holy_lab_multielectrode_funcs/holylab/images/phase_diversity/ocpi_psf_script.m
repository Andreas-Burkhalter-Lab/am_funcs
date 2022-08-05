%------------------------------------------
imsz = [256 256]; % Size of image
pixel_spacing = [0.1 0.1]; %microns

zdim = 30; % the thickness of stack, inmicrons
zspace = 0.1; % z-pixel spacing, in microns

zFrame = -(zdim)/2:zspace:(zdim)/2; % z-depth vector

%------------------------------------------
Zcoefs = zeros(1,9); % zernike aberration vector, starts from defocus, Zcoef(1) = defocus   
%Zcoefs(4) = -6;

%-------------------------------------------

n = 1.33; % refractive index of water
lambda_det = 0.53; % wavelength of light in free space
NA_det = 0.5; % NA of detector (objective lens)

% The illumination parameters are tuned for OCPI.2
lambda_ill = 0.488;
NA_grin = 0.2; % NA of full cylindrical lens
laser_e2 = 700.0; % laser beam diameter, in microns
slit_x = 300; % slit width at back aperture of grin lens, in microns
grin_thickness = 1000; % total thickness of grin lens, in microns

%----------------------------------------------
%----------------------------------------------



[psf_fm, psf_cyl, psf_ocpi] = ocpi_psf(imsz, zFrame, pixel_spacing, Zcoefs, n, lambda_det, NA_det, lambda_ill, NA_grin, laser_e2, slit_x, grin_thickness);

% psf_fm_temp = permute(psf_fm, [3 1 2]);
% psf_tot = sum(sum(psf_fm_temp(:,:,:),3),2);
% 
% figure; plot(psf_tot)


psf_fm = psf_fm/max(psf_fm(:));
psf_cyl = psf_cyl/max(psf_cyl(:));
psf_ocpi = psf_ocpi/max(psf_ocpi(:));
% 
% psf_fm = psf_fm.^0.4;
% psf_ocpi = psf_ocpi.^0.4;
% 
clims = [0 1];

%plotting the results in a useful way
figure;
mid_x = round(imsz(1)/2);
mid_y = round(imsz(2)/2);
nFrame = length(zFrame);
mid_z = round(nFrame/2);


psf_cyl_3d = repmat(psf_cyl, [1 1 size(psf_ocpi,3)]);
subplot(2,2,1); imshow((squeeze(psf_cyl_3d(mid_x, :,:))), clims);
psf_cyl_1d = psf_cyl(:, mid_x);
psf_cyl_xaxis = (1:nFrame)*zspace;
psf_cyl_fwhm = fwhm(psf_cyl_xaxis, psf_cyl_1d);
title(['Illumination: FWHM = ' num2str(psf_cyl_fwhm)]);

subplot(2,2,2); imshow((squeeze(psf_fm(mid_x, :,:)))', clims);
psf_fm_axial = squeeze(squeeze(psf_fm(mid_x, mid_y, :)));
psf_axial_xaxis = (1:nFrame)*zspace;
psf_axial_fwhm = fwhm(psf_axial_xaxis,psf_fm_axial);
title(['Detector: Axial FWHM = ' num2str(psf_axial_fwhm)]);


subplot(2,2,3); imshow(psf_ocpi(:,:,mid_z), clims);
psf_ocpi_lat = squeeze(psf_ocpi(:,mid_y,mid_z));
psf_lat_xaxis = (1:imsz(1))*pixel_spacing(1);
psf_lat_fwhm = fwhm(psf_lat_xaxis, psf_ocpi_lat);
title(['OCPI: Lateral FWHM = ' num2str(psf_lat_fwhm)]);


subplot(2,2,4); imshow((squeeze(psf_ocpi(mid_x, :,:)))', clims);
psf_ocpi_axial = squeeze(squeeze(psf_ocpi(mid_x, mid_y, :)));
psf_axial_xaxis = (1:nFrame)*zspace;
psf_axial_fwhm = fwhm(psf_axial_xaxis,psf_ocpi_axial);
title(['OCPI: Axial FWHM = ' num2str(psf_axial_fwhm)]);






