function [t,ro] = opt2dparabolicmirror(params,r)
% OPT2DMIRROR: defines a parabolic mirror for ray tracing
% very crude...DONT USE
% Syntax:
%   hline = opt2dmirror(params)
%   [t,ro] = opt2dmirror(params,ray)
% The first syntax plots the "surface" on the current axis. The second
% traces a ray.
%

  if (nargin == 1)
    % Plot surface on screen  
    y = linspace(-params.ymax,params.ymax,100);
    x = -params.curv * y.^2;
    t = line(x,y,'Color','k');
  else
    ro = r;           
    for i = 1:length(r)
        % This assumes the rays are horizontal!
        y = r(i).x0(2);
        x = -params.curv * y^2;
        t(i) = x-r(i).x0(1);
        l = sqrt(1+4*params.curv^2*y^2);
        sth = 1/l;
        cth = -2*params.curv*y/l;
        ro(i).x0 = [-params.curv*y^2 y];
        ro(i).e = [cth^2-sth^2, 2*sth*cth];
    end
  end
  % mirror