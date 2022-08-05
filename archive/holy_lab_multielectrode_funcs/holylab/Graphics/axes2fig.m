function cor_fig=axes2fig(cor_axes,hAxes)
% axes2fig: translate coordinates in axes to coordinates in figure
% syntax: 
%    cor_fig=axes2fig(cor_axes,hAxes)
% pre:
%    cor_axes: the coordinates(x,y) in data units;
%    hAxes=gca: the axes' handle
% post:
%    cor_fig: the coordinates in figure's units
% see: fig2axes, data2norm
% note: this is a wrapper to data2norm.
   if(nargin==1) hAxes=gca; end
   [cor_fig(1), cor_fig(2)]=data2norm(cor_axes(1), cor_axes(2),hAxes);
   pos=get(get(hAxes, 'parent'), 'Position');
   cor_fig=cor_fig.*pos(3:4);
   
