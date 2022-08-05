function [nOut,xCenterOut,yCenterOut] = hist2d(xIn,yIn,rect,nx,ny)
% HIST2D: 2-dimensional histograms
% Syntax:
%   hist2d(xIn,yIn,rect,nx,ny)
% plots the histogram as a grayscale intensity plot with log scaling
%   [nOut,xCenter,yCenter] = hist2d(xIn,yIn,rect,nx,ny)
% instead returns the numbers per bin.
%
% Bin data points given by (xIn,yIn)
% into bins with centered on xCenters,yCenters
% rect determines the exterior range:
% rect = [xmin xmax ymin ymax]
if (length(xIn) ~= length(yIn))
        error('x & y must have the same length!');
end
xBoundary = linspace(rect(1),rect(2),nx+1);
yBoundary = linspace(rect(3),rect(4),ny+1);
xCenter = (xBoundary(1:nx)+xBoundary(2:nx+1))/2;
yCenter = (yBoundary(1:ny)+yBoundary(2:ny+1))/2;
wx = (rect(2)-rect(1))/nx;
wy = (rect(4)-rect(3))/ny;
xi = round((xIn-xCenter(1))/wx)+1;
xOK = find(xi <= nx & xi >= 1);
yi = round((yIn-yCenter(1))/wy)+1;
yOK = find(yi <= ny & yi >= 1);
indx = intersect(xOK,yOK);
n = zeros(nx,ny);
for i = 1:length(indx)
        n(xi(indx(i)),yi(indx(i))) = n(xi(indx(i)),yi(indx(i)))+1;
end
if (nargout == 0)
        imagesc(xCenter,yCenter,log(n+1)');
        set(gca,'YDir','normal');
        colormap(1-gray);
else
        nOut = n;
        xCenterOut = xCenter;
        yCenterOut = yCenter;
end
