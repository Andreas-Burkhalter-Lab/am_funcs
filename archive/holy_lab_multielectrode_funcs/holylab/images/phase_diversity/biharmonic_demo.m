th = 2*pi*rand;
s = 2*[1 rand];
R = 6;
xi0 = R*(rand(1,2)-0.5);
m = 1;
zeta = 0;
NA = 0.5;
lambda = 0.5;
M = 20;
pixel_spacing = 4;
alpha = 0;  % defocus

im = double(imread('cameraman.tif'));

rot = [cos(th), -sin(th); sin(th) cos(th)];
A = rot*diag(s)*rot';
a = A([1 2 4]);
A = [1 0.2; -0.1 1.1];
a = A(:);
[act_grid,act_xi] = mirao52_layout;
[H0,rho,theta] = pupil_initialize(NA,M,lambda,size(im),pixel_spacing);
pupildata = struct('H0',H0,'rho',rho,'theta',theta);

param = struct('m',m,'zeta',zeta,'a',a,'xi0',xi0,'R',R,'act_xi',act_xi,'alpha',alpha);

v = zeros(2,52);
%v(2,35) = 40;
%v(2,30) = -40;
v(2,:) = 8*randn(1,52);

phik = biharmonic_penalty(param,v,pupildata);
[Hk,sk,im_obs] = pd_forward_model_2d(phik,H0,im);

%param0 =
%struct('m',1.1*m,'zeta',zeta,'a',a*0.9,'xi0',xi0*1.05,'R',R*1.3,'act_xi',act_xi,'alpha',1.2*alpha);
param0 = param;
param0.a = param0.a + 0.02*randn(size(param0.a));

% Set up optimization of particular parameters
%extract_func = @(s) extract_fields(s,'m','zeta','a','xi0','R','alpha');
extract_func = @(s) extract_fields(s,'a');
[p0,fields,field_shape,sbase] = extract_func(param0);
fill_func = @(p) fill_fields(fields,field_shape,p,sbase);
myfunc = @(s) biharmonic_penalty(s,v,pupildata,im_obs);
opt_func = @(p) optimize_struct_wrapper(p,myfunc,extract_func,fill_func);

% Do the optimization
p = fminunc(opt_func,p0,optimset('GradObj','on','Display','iter'));
param1 = fill_func(p);
