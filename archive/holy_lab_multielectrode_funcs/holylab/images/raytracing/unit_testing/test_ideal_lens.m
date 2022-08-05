function test_ideal_lens
lens.nL = 1.5;
lens.nR = 1;
lens.zL = 5;
lens.zR = 10;
lens.fL = -5;
lens.fR = 20;

n_rays = 21;
thmax = 0.5;
ymax = 0.9;

th = linspace(-1,1,n_rays)*thmax;
y = linspace(-1,1,n_rays)*ymax;

% An on-axis, in-focal-plane point
x0 = zeros(3,n_rays);
s0 = [0*th; sin(th); cos(th)];
% [x1,s1] = raytrace_ideal_lens(x0,s0,lens);
% il_plot(x0,s0,1,x1,s1,1);
% 
% % Trace in reverse direction
% [x0r,s0r] = raytrace_ideal_lens(x1,-s1,lens);
% il_plot(x1,s1,1,x0r,s0r,1);
% 
% s0+s0r  % this should be zero within roundoff
% 
% % An off-axis point
% x0 = repmat([0; 0.2; 0],1,n_rays);
% [x1,s1] = raytrace_ideal_lens(x0,s0,lens);
% il_plot(x0,s0,1,x1,s1,1);
% % Reverse it
% [x0r,s0r] = raytrace_ideal_lens(x1,-s1,lens);
% il_plot(x1,s1,1,x0r,s0r,1);
% s0+s0r
% 
% % Trace from out-of-focus origins (but rays still converge to a point in
% % the object focal plane)
% x0a = x0 - s0/2;
% [x1,s1] = raytrace_ideal_lens(x0a,s0,lens);
% il_plot(x0a,s0,1,x1,s1,1);
% 
% % Infinity-corrected lenses, simple model
% lensinf = lens;
% lensinf.fR = inf;
% [x1,s1] = raytrace_ideal_lens(x0,s0,lensinf);
% il_plot(x0,s0,1,x1,s1,1);
% % Reverse it
% [x0r,s0r] = raytrace_ideal_lens(x1,-s1,lensinf);
% il_plot(x1,s1,1,x0r,s0r,1);
% s0+s0r
% 
% % Other side
% lensinf = lens;
% lensinf.fL = -inf;
% x0 = zeros(3,n_rays); x0(2,:) = y;
% s0 = zeros(3,n_rays); s0(2,:) = 0.1;
% s0(3,:) = sqrt(1-sum(s0(1:2,:).^2,1));
% [x1,s1] = raytrace_ideal_lens(x0,s0,lensinf);
% il_plot(x0,s0,1,x1,s1,1);
% % Reverse it
% [x0r,s0r] = raytrace_ideal_lens(x1,-s1,lensinf);
% il_plot(x1,s1,1,x0r,s0r,1);
% s0+s0r
% 
% 
% % Have the general model mimic the simple one by supplying an explicit z_N
% % that gives ell = f.
% phi = linspace(0,0.05,21);
% tmp = 1-cos(phi);
% nodes = [phi; lens.fL*tmp; lens.fR*tmp];
% lensg = prepare_ideal_lens([lens.fL lens.fR],nodes,[],[]);
% lensg = copyfields(lens,{'nL','nR','zL','zR'},lensg);
% x0 = repmat([0; 0.2; 0],1,n_rays);
% s0 = [0*th; sin(th); cos(th)];
% [x1,s1] = raytrace_ideal_lens(x0,s0,lens);
% [x1g,s1g] = raytrace_ideal_lens(x0,s0,lensg);
% [x0g,s0g] = raytrace_ideal_lens(x1,-s1,lensg);
% s0g+s0
% 
% % Note: this next one is hard to mimic in the general model. It requires
% % that z_N(phi)/ell -> 1-cos(phi) as ell->inf. Not really feasible with the
% % infrastructure developed here.
% % % An infinite-conjugate by the general model
% % lensinf = lens;
% % lensinf.fR = inf;
% % nodes = [phi; lensinf.fL*tmp; 0*phi]; % here's the oddity...
% % lensg = prepare_ideal_lens([lensinf.fL lensinf.fR],nodes,[],[]);
% % lensg = copyfields(lensinf,{'nL','nR','zL','zR'},lensg);
% % [x1,s1] = raytrace_ideal_lens(x0,s0,lensinf);
% % [x1g,s1g] = raytrace_ideal_lens(x0,s0,lensg);
% % s1-s1g
% % x1-x1g

% Do the aplanatic sphere
n_in = 1.6;
n_out = 1.2;
nratio = n_in/n_out;
R = 1;
R_in = R/nratio;
R_out = R*nratio;
ang = linspace(0,1,51)*(0.99*pi/2);
fc = [sin(ang); 1-cos(ang)];
lens_aplanat = prepare_ideal_lens(-[R_in R_out],[],R_in*fc,R_out*fc);
lens_aplanat.nL = n_in;
lens_aplanat.nR = n_out;
lens_aplanat.zL = 0;
lens_aplanat.zR = 0;
%ang0 = -0.1;
ang0 = 0;
x0 = repmat([0; -sin(ang0); -cos(ang0)]*R_in, 1, n_rays);
[x1,s1] = raytrace_ideal_lens(x0,s0,lens_aplanat);
il_plot(x0,s0,2,x1,s1,2.5);
% Calculate the place of intersection
for i = 1:n_rays
  A = [s0(2:3,i) -s1(2:3,i)];
  b = x1(2:3,i) - x0(2:3,i);
  t = A\b;
  xc(:,i) = x0(2:3,i) + t(1)*s0(2:3,i);
  xc2(:,i) = x1(2:3,i) + t(2)*s1(2:3,i);
end
rad2 = sum(xc.^2,1);
rad2  % should all be at radius 1
% Try backwards
[x0a,s0a] = raytrace_ideal_lens(x1,-s1,lens_aplanat);
il_plot(x1,s1,2.5,x0a,-s0a,2);
% Try it with rays not starting on surface
x0n = x0-0.5*s0;
[x1n,s1n] = raytrace_ideal_lens(x0n,s0,lens_aplanat);
il_plot(x0n,s0,2.5,x1,s1,2.5);


% Try it with the lens flipped
fc(2,:) = -fc(2,:);
lens_aplanat = prepare_ideal_lens([R_out R_in],[],R_out*fc,R_in*fc);
lens_aplanat.nL = n_out;
lens_aplanat.nR = n_in;
lens_aplanat.zL = 0;
lens_aplanat.zR = 0;
ang0 = -0.2;
x0 = repmat([0; -sin(ang0); cos(ang0)]*R_in, 1, n_rays);
s0 = -s0;
[x1,s1] = raytrace_ideal_lens(x0,s0,lens_aplanat);
il_plot(x0,s0,2,x1,s1,2.5);
% Try backwards
[x0a,s0a] = raytrace_ideal_lens(x1,-s1,lens_aplanat);
il_plot(x1,s1,2.5,x0a,-s0a,2);

% Test inffoc: model a Luneburg lens
n_out = 1;
R = 1;
ang = linspace(0,1,51)*(0.99*pi/2);
fc = [sin(ang); 1-cos(ang)];
lens_lune = prepare_ideal_lens([-inf R],[],[],R*fc);
lens_lune.nL = n_out;
lens_lune.nR = n_out;
lens_lune.zL = 0;
lens_lune.zR = 0;
x0 = zeros(3,n_rays); x0(2,:) = y; x0(3,:) = -2;
s0 = zeros(3,n_rays); s0(2,:) = 0; s0(3,:) = sqrt(1-sum(s0(1:2,:).^2));
[x1,s1] = raytrace_ideal_lens(x0,s0,lens_lune);
il_plot(x0,s0,2,x1,s1,2.5);
% Trace backwards
[x0a,s0a] = raytrace_ideal_lens(x1,-s1,lens_lune);
il_plot(x0,s0,2,x1,s1,2.5);

function hline = il_plot(varargin)
  figure
  for i = 0:3:3
    x = varargin{1+i}([3 2],:);
    s = varargin{2+i}([3 2],:);
    t = varargin{3+i};
    hline(i) = line([x(1,:); x(1,:)+t*s(1,:)], [x(2,:); x(2,:)+t*s(2,:)]);
  end
  