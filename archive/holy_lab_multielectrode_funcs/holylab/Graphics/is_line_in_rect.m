function result=is_line_in_rect(hLine, rect)
   % NOTE: maybe a better one is is_line_cross_rect.
   % NOTE: rect is in the form of (x,y, width, height)
   % SEE: is_in_rect()
   
   x=get(hLine, 'xdata');
   y=get(hLine, 'ydata');
   if(any(x>=rect(1) & y>=rect(2) & x<=rect(1)+rect(3) & y<=rect(2)+rect(4)))
      result=1;
   else
      result=0;
   end
   
