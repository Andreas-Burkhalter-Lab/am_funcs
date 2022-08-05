function plotSeparationLines(hAxes, separationPositions)
% plot vertical dotted lines

   ylim=get(hAxes, 'ylim');
   hold(hAxes, 'on');
   if(~isempty(separationPositions))
      hDotLines=plot(repmat(separationPositions, 2,1), ylim');
      set(hDotLines, 'lineStyle', ':', 'color', 'black', 'hitTest', 'off');
   end
   
