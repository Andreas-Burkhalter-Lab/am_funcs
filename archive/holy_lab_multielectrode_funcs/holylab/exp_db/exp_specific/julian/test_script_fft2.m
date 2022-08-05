%(playing around script) trying to get my head around fft2 and masked
%images

entry = load('/usr/lab/exp_db/julian/julian_141886.xdb', '-mat');

[masked_imgs, masks] = jpm_region_mask(entry);

fftd_img(1) = {fft2(masked_imgs{1}{1})};
fftd_img(2) = {fft2(masked_imgs{1}{2})};

fftd_mask(1) = {fft2(masks{1}{1})};
fftd_mask(2) = {fft2(masks{1}{2})};

ratiod_img(1) = {fftd_img{1}./fftd_mask{1}};
ratiod_img(2) = {fftd_img{2}./fftd_mask{2}};

subd_img(1) = {fftd_img{1}-fftd_mask{1}};
subd_img(2) = {fftd_img{2}-fftd_mask{2}};

radial_dists = [];

middle_x_pixel = floor(size(fftd_img{1},2)/2);
middle_y_pixel = floor(size(fftd_img{1},1)/2);

figure(3); hold on;
for idx = 1:size(subd_img{1},1)
    for idy = 1:size(subd_img{1},2)
        radial_dists(idx,idy) = sqrt((idx-middle_x_pixel)^2+(idy-middle_y_pixel)^2);
    
    end
    plot(radial_dists(idx,:),abs(real(subd_img{1}(idx,1:size(radial_dists(idx,:),2)))));
end
set(gca,'yscale', 'log');

figure(4); hold on;
for idx = 1:size(subd_img{2},1)
    for idy = 1:size(subd_img{2},2)
        radial_dists(idx,idy) = sqrt((idx-middle_x_pixel)^2+(idy-middle_y_pixel)^2);
    
    end
    plot(radial_dists(idx,:),abs(real(subd_img{2}(idx,1:size(radial_dists(idx,:),2)))));
end
set(gca,'yscale', 'log');