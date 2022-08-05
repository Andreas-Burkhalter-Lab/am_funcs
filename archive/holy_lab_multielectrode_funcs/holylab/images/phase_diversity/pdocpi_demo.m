clear sigma_history Zc_history
imclass = 'double';  % double good for numerically testing gradients
plotflag = true;
gradtest = true;

% Use the next 3 lines for a "real" image
immid = cast(imread('cameraman.tif'),imclass);
imsz = [size(immid) 10]; mid = round(imsz(end)/2);
im = zeros(imsz,imclass);
im(:,:,mid) = immid;
% Or, uncomment the next lines for a single "bead"
imsz = [32 32 40];
im = zeros(imsz,imclass);
for indx = 1:length(imsz)
 imc{indx} = round(imsz(indx)/2);
end
mid = imc{end};
im(imc{:}) = 1;

%noisefactor = 0.01;
noisefactor = 0;

% Set up the imaging conditions
NAred = 0.5/1.34;
M = 20;
lambda = 0.5;
defoc = [0 20*lambda];
K = length(defoc);
pixel_spacing = [4 4];
zstep = 0.5;  % spacing between frames in stack in microns
zw = 1.5;  % sigma of illumination in microns
% stacklen = imsz(3);  % # of images in stack
% z = 0:stacklen-1;
% zwrap = (z > stacklen/2);
% z(zwrap) = z(zwrap) - stacklen;
% v = exp(-(z*zstep/zw).^2/2);
% v = v / sum(v);

[H0,rho,theta,c0] = pupil_initialize(NAred,M,lambda,size(im),pixel_spacing,zstep);
H0 = cast(H0,imclass);

% Define a coma aberration and the amount of defocus
Z7coef = 3;
Zinfo = struct('Zindex',{{7,[4 7]}},'Zcoefs',{{Z7coef,[defoc(2) Z7coef]}},'map',{{1, [0 1]}});
% Set up the structure for computation
pdparams = struct('rho',rho,'theta',theta,'H0',H0,'c0',c0,'sigma',zw/zstep,...
  'Zinfo',Zinfo,'f',im);
outputparams = struct('PSFs',false);
% Compute the observed images
imdata = ocpi_pdwrapper(pdparams,outputparams);
pdparams.Iobsk = imdata.Iobs_calc; % Store them as targets for minimization

% Now insert some wrong parameters, so that we have something to optimize
pdparams.f = imdata.Iobs_calc(:,:,:,1);
pdparams.sigma = 1;
for indx = 1:K
  map = pdparams.Zinfo.map{indx};
  Zindx = find(map > 0);  % zero anything that will be optimized
  pdparams.Zinfo.Zcoefs{indx}(Zindx) = 0;
end

if plotflag
  % x-y plots
  n_from_mid = 2;
  skip = 2;
  rng = mid-n_from_mid*skip:skip:mid+n_from_mid*skip;
  plot_stack_row(im(:,:,rng));
  suptitle('Original object')
  plot_stack_row(imdata.Iobs_calc(:,:,rng,1));
  suptitle('Camera 1')
  plot_stack_row(imdata.Iobs_calc(:,:,rng,2));
  suptitle('Camera 2')
end

if gradtest
  s = ocpi_pdwrapper(pdparams,struct('fgrad',true,'vgrad',true,'phikgrad',true));
  % Check one of the f derivatives
  midfIndex = round(numel(pdparams.f)/2);
  df = 1e-8;
  pdparams1 = pdparams;
  pdparams1.f(midfIndex) = pdparams.f(midfIndex) + df;
  s1 = ocpi_pdwrapper(pdparams1);
  fprintf('f grad: analytic %g, numeric %g\n',s.fgrad(midfIndex),(s1.val-s.val)/df);
  % Check the sigma derivative
  dsigma = 1e-8;
  pdparams1 = pdparams;
  pdparams1.sigma = pdparams.sigma + dsigma;
  s1 = ocpi_pdwrapper(pdparams1);
  fprintf('sigma grad: analytic %g, numeric %g\n',s.sigmagrad,(s1.val-s.val)/dsigma);
end

% Let's pre-optimize the object, so that it doesn't wildly throw off the
% other parameters, by taking 10 conjugate-gradient steps
fprintf('Pre-optimizing the object\n');
pdparams.datafields = {'f'};
myfun = @(p) ocpi_pdwrapper(p,pdparams);
p = ocpi_pdwrapper(pdparams,'extract');
[p,err0] = conjgrad(myfun,p,struct('Display',true,'iter_max',10));
pdparams = ocpi_pdwrapper(p,pdparams,'insert');

% Optimize all the parameters
pdparams.datafields = {'f','sigma','Zcoefs'};
p = ocpi_pdwrapper(pdparams,'extract');
err = err0(end);
mu = ones(1,length(pdparams.datafields));
sigma_history = pdparams.sigma;
Zc = ocpi_pdwrapper(pdparams,'extract',{'Zcoefs'})'; % for history
Zc_history = Zc;
if plotflag
  % We'll have a continuous view of the progress in terms of parameters
  figure; hline_sigma = plot(sigma_history); title('sigma')
  figure; hline_Zc = plot(Zc_history); title('Zernikes')
end
for i = 1:100
  [pdparams,err(i),mu] = pdimprove_all(pdparams,mu);
  sigma_history(end+1) = pdparams.sigma;
  Zc = ocpi_pdwrapper(pdparams,'extract',{'Zcoefs'});
  Zc_history(end+1,:) = Zc;
  if plotflag
    set(hline_sigma,'YData',sigma_history);
    set(hline_Zc,'YData',Zc_history);
    drawnow
  end
end

if plotflag
  plot_stack_row(pdparams.f(:,:,rng));
  suptitle('Final object')
end

return

%p(end+1:end+3) = [0; defoc(2); 0];
%pdparams.datafields = {'f','sigma','Zcoefs'};
myfun = @(p) ocpi_pdwrapper(p,pdparams);

[p,err] = conjgrad(myfun,p,struct('Display',true,'iter_max',500,'goodstep_fcn',@(p) phistory(p,length(p),'sigma_history')));

pdopt = ocpi_pdwrapper(p,pdparams,'insert');

if plotflag
  figure
  semilogy(err0,'r'); hold on; semilogy(err,'b'); hold off
  figure; plot(sigma_history); title('sigma')
  plot_stack_row(pdopt.f(:,:,rng))
end

return

Zdefoc = zernike_values(rho,theta,4);
% Define a phase aberration
Zval = zernike_values(rho,theta,7); % a Z(3,-1) aberration
phi = 1*Zval;
phik = zeros([imsz(1:2) K],imclass);
for Kindx = 1:K
  phik(:,:,Kindx) = phi + defoc(Kindx) * Zdefoc;
end

% Generate the observed images
s = phasedivocpi_val_grad(phik,v,H,rho,c0,sqrt(im));
Iobsk = s.Iobs;
mask = ones(size(im),imclass);


if gradtest
  if ~exist('psi0','var')
    psi0 = cast(rand(imsz),imclass);
  end
  psi0 = psi0/sqrt(sum(psi0(:).^2)/sum(im(:))); % Give the same total brightness
  phik0 = zeros([imsz(1:2) K],imclass);
  opst = struct('psigrad',true,'vgrad',true,'phikgrad',true);
  opsf = struct('psigrad',false,'vgrad',false,'phikgrad',false);
  s = phasedivocpi_val_grad(phik0,v,H,rho,c0,psi0,Iobsk,mask,opst);
  
  return
  % Check the psi gradient
  dpsi = 1e-8;
  psigrad = s.psigrad;
  for c1 = 1:imsz(1)
    for c2 = 1:imsz(2)
      for c3 = 1:imsz(3)
        psi1 = psi0;
        psi1(c1,c2,c3) = psi0(c1,c2,c3) + dpsi;
        stest = phasedivocpi_val_grad(phik0,v,H,rho,c0,psi1,Iobsk,mask,opsf);
        psigrad(c1,c2,c3) = (stest.val - s.val)/dpsi;
      end
    end
  end
  
  % Check the v gradient
  dv = 1e-8;
  vgrad = s.vgrad;
  for c1 = 1:length(v)
    v1 = v;
    v1(c1) = v(c1) + dv;
    stest = phasedivocpi_val_grad(phik0,v1,H,rho,c0,psi0,Iobsk,mask,opsf);
    vgrad(c1) = (stest.val - s.val)/dv;
  end
  
  % Check the phi gradient
  dphi = 1e-8;
  phikgrad = s.phikgrad;
  for c1 = 1:size(phikgrad,1)
    for c2 = 1:size(phikgrad,2)
      for c3 = 1:size(phikgrad,3)
        phik1 = phik0;
        phik1(c1,c2,c3) = phik0(c1,c2,c3) + dphi;
        stest = phasedivocpi_val_grad(phik1,v,H,rho,c0,psi0,Iobsk,mask,opsf);
        phikgrad(c1,c2,c3) = (stest.val - s.val)/dphi;
      end
    end
  end
end
