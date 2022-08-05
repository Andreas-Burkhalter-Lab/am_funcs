% Use this next line for a "real" image
im = double(generate_fake_3d_stack);
im = imreduce(single(im), [4 4 1]);
% im = double(im);
%   im = zeros([128 128 128]);
%   im(64, 64, 64) = 1;

% im = zeros([64 64 64]);
% im(32, 32, 32) = 256;
% fname = '/mnt/dturaga_005/dturaga/phase_diversity_3d/tyr_MPR_native.4dfp.img';
% fid = fopen(fname,'r','ieee-be');
% img = double((fread(fid,'float32')));
% fclose(fid);
% im = reshape(img,144,172,120);
% %im = imreduce(single(im), [2 2 2]);
% im = double(im);
% 
% 
% noisefactor = 0.0005;
noisefactor = 0.00;

%------------------------------------------
imsz = size(im); % Size of image
pixel_spacing = [0.25 0.25 0.5]; %microns
%pixel_spacing = [0.2 0.2 0.7];

zFrame = -(imsz(3)*pixel_spacing(3))/2:pixel_spacing(3):((imsz(3)-1)*pixel_spacing(3))/2; % z-depth vector


%-------------------------------------------

n = 1.33; % refractive index of water
lambda_det = 0.53; % wavelength of light in free space
NA_det = 0.5; % NA of detector (objective lens)
%NA_det = 0.8;


defoc = 7*lambda_det; % magnitude of known aberration
[H,pdk,rho,theta] = pdinitialize3d(NA_det, lambda_det, n, imsz, pixel_spacing,[0 defoc]);

%----------------------------------------------
% calculating 3d psf and otf
% adapted from hanser et al, J. Microscopy, vol 216, pg 32, year 2004.

coef = n/lambda_det;
kz = (coef^2 - rho.^2).^0.5;

zFrame = abs(zFrame);

Hz = zeros(imsz, class(H));
hz = zeros(imsz, class(H));

for indx = 1:imsz(3)
    Hz(:,:,indx) = H .* exp(i * 2 * pi * zFrame(indx) .* kz);
    hz(:,:,indx) = ifft2(Hz(:,:,indx)); 
end


h2 = hz .* conj(hz); % 3d psf
H2 = fftn(h2);
imfft = fftn(im);

imdl = ifftn(imfft .* H2);

% 
% figure;
% % for indx = 1:round(imsz(3)/10):imsz(3)
% for indx = 1:imsz(3)
%     suptitle(num2str(indx));
%     subplot(2,2,1); imshowsc(im(:,:,indx)); axis image; title('Object')
%     subplot(2,2,2); imshowsc(imdl(:,:,indx)); axis image; title('Diffraction-limited image')
%     subplot(2,2,3); imshowsc(squeeze(im(indx, :,:))); title('Object-z');
%     subplot(2,2,4); imshowsc(squeeze(imdl(indx, :,:))); title('Object-z');
%     pause;
%     
% end
% 
% temp_imdl = permute(imdl, [3 1 2]);
% temp_imdl = sum(sum(temp_imdl(:,:,:),3),2);
% figure; plot(temp_imdl);


%------------------------------------------------------------------------
% Calculating the two phase aberrated images

% Define a phase aberration
rholim = rho;
rholim(rho > 1) = 0;
phi = reshape(zernfun2(7,rholim(:),theta(:)),size(rho)) .* H; % a Z(3,-1) aberration
phi = double(phi);
phi = 3 * phi;
%phi = 0 * phi;

zFrame = abs(zFrame);

%figure; imagesc(phi); title('Phase aberration'); colorbar

imk = zeros([size(Hz), 2],class(im));
K = size(pdk,3);


Hk = zeros(size(H), class(H));


Hzk = zeros(size(Hz), class(H));
hzk = zeros(size(Hz), class(H));


for indx1 = 1:K
    Hk(:,:,indx1) = H .* exp(i * (phi + pdk(:,:,indx1)));
    
    for indx2 = 1:imsz(3)
        
        Hzk(:,:,indx2, indx1) = Hk(:,:,indx1) .* exp(i * 2 * pi * zFrame(indx2) .* kz);
        hzk(:,:,indx2) = ifft2(Hzk(:,:,indx2, indx1));   
        
    end
    
    hzk2 = hzk .* conj(hzk); 
    Hzk2 = fftn(hzk2);
    imfft = fftn(im);
    imk_temp = ifftn(imfft .* Hzk2);
    imk(:,:,:, indx1) = imk_temp;
end

%---------------------------------------
%adding noise

if noisefactor ~= 0
    for indx = 1:size(imk,4)
        tmp = imk(:,:,:,indx);
        noiseamplitude = noisefactor * mean(tmp(:));
        imk(:,:,:,indx) = imk(:,:,:,indx) + noiseamplitude*randn(size(tmp));
    end
end


% figure;
% for indx = 1:imsz(3)
%     
%     suptitle(num2str(indx));
%     
%     subplot(3,3,1); imshowsc(im(:,:,indx)); axis image; title('original image');
%     subplot(3,3,4); imshowsc(squeeze(im(indx,:,:))); title('original image -z');
%     
%     subplot(3,3,2); imshowsc(imdl(:,:,indx)); axis image; title('diffraction limited');
%     subplot(3,3,3); imshowsc(squeeze(imdl(indx, :,:))); title('diffraction limited - z');
%     
%     subplot(3,3,5); imshowsc(imk(:,:,indx,1)); axis image; title('abr 1');
%     subplot(3,3,8); imshowsc(imk(:,:,indx,2)); axis image; title('abr 2');
%     
%     subplot(3,3,6); imshowsc(squeeze(imk(indx, :,:,1))); title('abr 1 - z');
%     subplot(3,3,9); imshowsc(squeeze(imk(indx, :,:,2))); title('abr 2 -z');
%     
%     pause;
%     
% end


Zcoefs = zeros(1,8); % from Z(2,0) up to Z(4,0)
%Zcoefs(4) = 1.3;
[Zcoefs,phi] = pdopt_zernike3d(imk,H,pdk,rho,theta,Zcoefs, Hz);
Zcoefs

[val,grad,object] = pdpenalty3d(phi,imk,H,pdk, Hz);
% figure; imshowsc(object); title('Estimate of original object')
% figure; imagesc(phi); title('Estimate of aberration')

 imdl = fftshift(imdl,3);
 imk = fftshift(imk,3);


%object = imfilter_gaussian(single(object), [3 3 3]);

figure;
for indx = 1:size(im,3)
    
    suptitle(num2str(indx));
    
    subplot(4,3,1); imshowsc(im(:,:,indx)); axis image; title('original image');
    subplot(4,3,2); imshowsc(imdl(:,:,indx)); axis image; title('diffraction limited');
    subplot(4,3,4); imshowsc(imk(:,:,indx,1)); axis image; title('aberration 1')
    subplot(4,3,5); imshowsc(imk(:,:,indx,2)); axis image; title('aberration 2')
    subplot(4,3,3); imshowsc(object(:,:,indx)); axis image; title('reconstructed')
    
    subplot(4,3,7); imshowsc(squeeze(im(indx,:,:))); axis image; title('original image');
    subplot(4,3,8); imshowsc(squeeze(imdl(indx,:,:))); axis image; title('diffraction limited');
    subplot(4,3,10); imshowsc(squeeze(imk(indx, :,:, 1))); axis image; title('aberration 1')
    subplot(4,3,11); imshowsc(squeeze(imk(indx,:,:,2))); axis image; title('aberration 2')
    subplot(4,3,9); imshowsc(squeeze(object(indx, :,:))); axis image; title('reconstructed')
    pause;
    
end

