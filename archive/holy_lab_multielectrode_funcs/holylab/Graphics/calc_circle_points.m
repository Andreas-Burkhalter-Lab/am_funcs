function [x,y]=calc_circle_points(def)
   roi.x=def(1); roi.y=def(2); roi.xyradius=def(3);
   npts = 100;
   th = linspace(0,2*pi,npts);
   x = roi.x + roi.xyradius*cos(th);
   y = roi.y + roi.xyradius*sin(th);
   
