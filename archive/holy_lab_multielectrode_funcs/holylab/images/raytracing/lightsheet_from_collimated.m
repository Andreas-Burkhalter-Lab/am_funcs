% A script for testing different lightsheet configurations, assuming we
% have as input a collimated beam (hopefully achromatized)

% Characteristics of the input beam
lambda = [440 488 515 561];
%lambda = [440];
diam = 1;
vertoffset = 0;
w = [1 Inf];  % width of aperture in || plane (limits horizontal field of view)

% Define the lens system
%ls = optical_stock_lenses('eo_32299'); csgap = 4.5;  % 12.5mm lens
ls = optical_stock_lenses('eo_32301'); csgap = 4.5; %csgap = 7;  % 15mm lens
%ls = optical_stock_lenses('eo_49925');  % 17.5mm lens
%lc = optical_stock_lenses('tl_LK1395L1');
lc = optical_stock_lenses('eo_48371');  % The 
contact_air = struct('t',0,'material','air');
cs = optical_stock_lenses('coverslip1.5');

% Layout for cylindrical-then-achromat
% layout = optical_layout_onaxis(contact_air,...
%   opticalAperture(w),...
%   struct('t',2,'material','air'),...
%   lc,...
%   contact_air,...
%   ls,...
%   struct('t',0,'material','water'));

% Layout for achromat-then-cylindrical
gap = lc{1}.span_along_optic_axis;  % determine how much spacing we need to avoid bumping (approximate)
gap = -gap(1);
layout_sequence = {contact_air,...
  opticalAperture(w),...
  struct('t',2,'material','air'),...
  ls,...
  struct('t',gap,'material','air'),...
  lc};
% Add a coverslip at the end, to decrease the working distance
throt = 0*pi/180;  % how much to rotate it to make it imperfectly aligned
qrot = [cos(throt/2), sin(throt/2)*[1 0 0]];  % rotation quaternion 
layout_sequence = [layout_sequence {...
  struct('t',csgap,'material','air'),...
  optical_rotate_component(cs,qrot)}];

% Immerse in water
layout_sequence = [layout_sequence {...
  struct('t',0,'material','water')}];

% Assemble the layout
layout = optical_layout_onaxis(layout_sequence{:});
nimmerse = opt_refrindx(layout{end},lambda);

% Run the tracing
wd = zeros(1,length(lambda));
dz = linspace(-5,0.2,1001)';  % for displacement from focal plane
negflag = dz < 0;
sigma2 = zeros(length(dz),length(lambda));
W0 = zeros(1,length(lambda));
z0 = zeros(1,length(lambda));
for i = 1:length(lambda)
  % Define the rays
  rb0 = optical_collimated_bundle([0 vertoffset 0],[0 0 1],lambda(i),3*diam,diam,[],51);
  % Trace all the way to focus along the [0 1 0] axis
  v = [0 1 0]';
  options = struct('tofocus',[[0 0 1]', v]);
  [rb,mean_coef,cov_coef] = optical_trace(rb0,layout,options);
  pf = mean_coef(:,2);
  c = [v'*cov_coef(:,:,1)*v,v'*cov_coef(:,:,2)*v,v'*cov_coef(:,:,3)*v];
  wd(i) = pf(3) - layout{end-1}.origin(3);
  % Analyze width as a function of displacement from focus
  sigma2(:,i) = c(1)*dz.^2 + c(2)*dz + c(3);
  % Calculate the illumination NA and gaussian beam parameters
  p = polyfit(dz(negflag),sqrt(sigma2(negflag,i)),1);
  theta_illum = atan(-p(1));
  NAillum = nimmerse(i)*sin(theta_illum);
  lambdamm = 1e-6*lambda(i); % lambda in millimeters
  W0(i) = lambdamm/(pi*nimmerse(i)*theta_illum);
  z0(i) = nimmerse(i)*pi*W0(i)^2/lambdamm;
  fprintf('lambda %g: working distance = %g, NAillum = %g, spot radius (no diffrac.) = %g, diffr. limited spot = %g, Rayleigh range %g\n',lambda(i),wd(i),NAillum,sqrt(c(3)),W0(i),z0(i))
end
fprintf('Working distance spread: ');
fprintf('%g ',wd-mean(wd));
fprintf('\n')

% Plot
figure
optical_plot_components_3d(layout)
optical_plot_rays_3d(rb)
figure
I = rb(end).intensity;
w = rb(end).wavelength;
col = squeeze(spectrumRGB(w(:)));
col = lightencolors(col,1-I(:));
p = rb(end).p;
scatter(p(1,:),p(2,:),16,col,'filled')

% Analyze intensity as a function of horizontal position in the sheet
clust = msams(p(1,:));
[clabel,nlabel] = agglabel(clust);
n_clust = length(clabel);
I = rb(end).intensity;
x = zeros(1,n_clust);
Ix = zeros(1,n_clust);
for i = 1:n_clust
  groupIndex = clabel{i};
  x(i) = mean(p(1,groupIndex));
  Ix(i) = sum(I(groupIndex));
end
[x,sortOrder] = sort(x);
Ix = Ix(sortOrder);
figure; plot(x,Ix); yl = get(gca,'YLim'); yl(1) = 0; set(gca,'YLim',yl)
xlabel('Horizontal position (mm)')
ylabel('Intensity (arb.)')

% Plot the width vs displacement
figure; hline = plot(dz,sqrt(sigma2));
for i = 1:length(hline)
  set(hline(i),'Color',spectrumRGB(lambda(i)));
end
hl = hline(1);
xlabel('dz'); ylabel('Sheet width'); title('Sheet width vs. focal displacement (no diffraction)');
% Plot width of diffraction-limited sheet (for reference)
W = bsxfun(@times,sqrt(1+bsxfun(@rdivide,dz,z0).^2),W0);
hline = line(dz,W,'LineStyle','--');
for i = 1:length(hline)
  set(hline(i),'Color',spectrumRGB(lambda(i)));
end
hl(2) = hline(1);
legend(hl,'Raytrace','Gaussian beam')
