obj = double(imread('cameraman.tif'));
imsz = size(obj);

NA = 0.5;
M = 20;
lambda = 0.5;
pixel_spacing = 6;

[H0,rho,theta] = pupil_initialize(NA,M,lambda,imsz,pixel_spacing);

phi = zeros([imsz 2],class(obj));

phi(:,:,2) = 5*H0.*rho.^2;

[Hk,sk,imk] = pd_forward_model_2d(phi,H0,obj);
figure; for indx = 1:2; subplot(1,2,indx); imshowsc(imk(:,:,indx)); end

[val,grad,obj] = pdpenalty(phi,imk,H0);

phi0 = H0.*(phi(:,:,2)+randn(size(rho)));
%phi0 = zeros(size(rho));
phiab = calcphi_split_path(imk,H0,phi0);

figure; subplot(1,3,1); imagesc(phi(:,:,2)); subplot(1,3,2); imagesc(phiab); subplot(1,3,3); imagesc(phi0)
