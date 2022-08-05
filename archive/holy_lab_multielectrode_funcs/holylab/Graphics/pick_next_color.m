function color=pick_next_color(aAxes)
   curIdx=getappdata(aAxes, 'colorIdx');
   if(isempty(curIdx))
      curIdx=1;
   else
      curIdx=curIdx+1;
   end
   colorOrder=get(aAxes, 'ColorOrder');
   curIdx=mod(curIdx, size(colorOrder,1));
   if(curIdx==0)
      curIdx=size(colorOrder,1);
   end
   color=colorOrder(curIdx, :);
   setappdata(aAxes, 'colorIdx', curIdx);
