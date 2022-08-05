function ip_result = imcroprect(ip,hCropRect)
   gd=guidata(hCropRect);
   x1=get(gd.line_x1, 'xdata');
   x2=get(gd.line_x2, 'xdata');
   y1=get(gd.line_y1, 'ydata');
   y2=get(gd.line_y2, 'ydata');
   
   ip_result=ip;
   for idx=1:length(ip_result)
      ip_result(idx).xrange=round(sort([x1(1),x2(1)]));
      ip_result(idx).yrange=round(sort([y1(1),y2(1)]));
   end
