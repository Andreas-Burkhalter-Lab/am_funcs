function cor_axes=fig2axes(cor_fig, hAxes)
% fig2axes: translate coordinates in figure to coordinates in axes
% syntax: 
%    cor_axes=fig2axes(cor_fig, hAxes)
% pre:
%    cor_fig: the coordinates in figure's units
%    hAxes=gca: the axes' handle
% post:
%    cor_axes: the coordinates(x,y) in data units;
% see: axes2fig
   
   if(nargin==1) hAxes=gca; end
   
   % convert cor_fig to be in unit "normalized"
   pos=get(get(hAxes, 'parent'), 'Position');
   cor_fig=cor_fig./pos(3:4);
   
   % Get position of axis in figure-normalized units
   axunits_old = get(hAxes,'Units');
   set(hAxes,'Units','normalized');
   axpos = get(hAxes,'Position');
   set(hAxes,'Units',axunits_old);
   
   % Get data limits
   xlim = get(hAxes,'XLim');
   ylim = get(hAxes,'YLim');
  
   % make this also works on log scale
   if(strcmp(get(gca, 'xscale'), 'log'))
      xlim=log(xlim);
   end
   if(strcmp(get(gca, 'yscale'), 'log'))
      ylim=log(ylim);
   end
  
   % Do the conversion
   cor_axes(1)=(cor_fig(1)- axpos(1))/axpos(3)*diff(xlim)+ xlim(1);
   cor_axes(2)=(cor_fig(2)- axpos(2))/axpos(4)*diff(ylim)+ ylim(1);
   
   if(strcmp(get(gca, 'xscale'), 'log'))
      cor_axes(1)=exp(cor_axes(1)); % 
   end
   if(strcmp(get(gca, 'yscale'), 'log'))
      cor_axes(2)=exp(cor_axes(2));
   end
