ng = [5 7];
z = round(20*rand(ng)+0.5);
im = zeros(30,100,70);
for i = 1:length(ng)
  xc{i} = round(linspace(1,size(im,1+i),ng(i)));
  xi{i} = 1:size(im,1+i);
end
Xc = cell(1,length(ng));
[Xc{:}] = ndgrid(xc{:});
Xi = cell(1,length(ng));
[Xi{:}] = ndgrid(xi{:});
zi = round(interpn(Xc{:},z,Xi{:}));
for i = 1:numel(zi)
  im(1:zi(i),i) = 1;
end

% Show the boundary
figure; surf(z); title('Real')
z0 = tissueboundary_initialize(im,0.5,[20 10]);
% How well are we doing so far?
figure; surf(z0); title('Initial guess')
figure; scatter(z(:),z0(:)); axis equal; title('Compare real vs. initial guess')
% Now try to improve
figure
p.height2intens = 0.5; imsz = size(im); p.gridsep = imsz(2:end)./ng; p.plot = true; p.maxiter = 1000;
p.z0 = z0;
[ztb,zitb] = tissueboundary(im,p);
title('Improved')
figure; scatter(z(:),ztb(:)); axis equal; title('Compare real vs. improved')
% Improve again
p.z0 = ztb;
figure
[ztb,zitb] = tissueboundary(im,p);
title('Improved2')
figure; scatter(z(:),ztb(:)); axis equal; title('Compare real vs. improved2')

