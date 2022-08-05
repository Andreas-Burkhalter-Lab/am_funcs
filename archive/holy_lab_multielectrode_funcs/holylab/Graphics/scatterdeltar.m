function hall = scatterdeltar(dr,drerr)
% SCATTERDELTAR: make 2-d deltar plots with errorbars
% Syntax:
%   scatterdeltar(dr,drerr)
%   linehandles = scatterdeltar(dr,drerr)
% where both inputs are 2-by-n matrices, with each row giving the change in
% firing rate (or, for drerr, the s.e.m. in firing rate change) to one
% particular stimulus and each column corresponding to a single cell or
% electrode.

% Copyright 2007 by Timothy E. Holy

  hy = errorbar(dr(1,:),dr(2,:),drerr(2,:));
  set(hy,'LineStyle','none');
  hold on
  hx = xerrorbar(dr(1,:),dr(2,:),drerr(1,:));
  delete(hx(2)); hx = hx(1);
  set([hx hy],'Color','k');
  hold off
  set(gca,'Box','off','TickDir','out')
  axis equal
  hall = [hx hy];
