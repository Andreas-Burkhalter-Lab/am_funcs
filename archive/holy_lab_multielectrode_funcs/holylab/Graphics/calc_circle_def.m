function def=calc_circle_def(hLine)
% calc_circle_def: calc the circle's definition given hLine's XData, YData
% pre: 
%    hLine: the handle of the circle
% post:
%    def: [x0, y0, r]
% see:
%    calc_circle_points()
   
   x=get(hLine, 'XData');
   y=get(hLine, 'YData');
   [x1,x2,x3]=deal(x(1), x(2), x(3));
   [y1,y2,y3]=deal(y(1), y(2), y(3));
   y0=calc_y0(x1,x2,x3,y1,y2,y3);
   x0=calc_y0(y1,y2,y3,x1,x2,x3); % note: x and y are symmetric
   r=((x1-x0)^2+(y1-y0)^2)^0.5;
   def=[x0 y0 r];
   
   
function y0=calc_y0(x1,x2,x3,y1,y2,y3)
   divided= ((x2*x2+y2*y2)-(x1*x1+y1*y1))*(x3-x1) ...
           -((x3*x3+y3*y3)-(x1*x1+y1*y1))*(x2-x1);
   
   divident= ( (x3-x1)*(y2-y1)-(x2-x1)*(y3-y1) )*2;
   
   y0=divided/divident;
   

