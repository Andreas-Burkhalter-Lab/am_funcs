function [t,ro] = opt2dDM(params,r)
% OPT2DMIRROR: defines a deformable mirror for ray tracing
% Syntax:
%   hline = opt2dmirror(params)
%   [t,ro] = opt2dmirror(params,ray)
% The first syntax plots the "surface" on the current axis. The second
% traces a ray.
%
% Params is the parameters structure describing the deformable surface:
%   mirEdge: edges of mirrors. 
% r is the input ray (see RAY)
% and
%   t is the distance along the ray to the point of intersection (NaN if no
%     intersection);
%   ro is the output ray after intersection.
%
% See also: RAY, OPT2DLINE, OPT2DCIRCLE, OPT2DMIRROR


  if (nargin == 1)
    % Plot surface on screen  
    t = line(params.mirEdge(:,1), params.mirEdge(:,2),'Color','k', 'LineWidth',2);
  else
    ro = r;           
    for j = 1:length(params.mirEdge)-1    
        % setting up defaults so as to be able to use the code from
        % opt2dmirror
        % check opt2dmirror for definitions
        x1 = params.mirEdge(j,1); y1 = params.mirEdge(j,2);
        x2 = params.mirEdge(j+1,1); y2 = params.mirEdge(j+1,2);
        params.mirPos = (x1+x2)/2;
        params.center = (y1+y2)/2;
        l = sqrt((y2-y1)^2 + (x2-x1)^2);
        sth = (y2-y1)/l;
        cth = (x2-x1)/l;
        params.length = l/2;
        % done
        params.normal = [-sth cth];     
        params.c = sum([params.mirPos params.center] .* params.normal);
        e_norm = sum(r.e.*params.normal);
        time = (params.c-sum(r.x0.*params.normal))/e_norm;
        p = r.x0 + time*r.e;  
        
        if (p(2) <= y2 && ...
              p(2) > y1 && ...
              time >= 0) % ie the ray hits this particular mirror
                
           ro.x0 = p;
           c2th = cth^2-sth^2;
           s2th = 2*sth*cth;
           ro.e = r.e *[ c2th s2th;...
                  s2th -c2th];
           t = time;
        end
     end
  end
  
  if ~exist('t') % ray didnt hit any of the mirrorlets
      t = NaN;
      ro = struct;
  end