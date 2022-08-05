function r = raytrace(r0,comp,drawflag)
% RAYTRACE: trace light rays through refractive surfaces
% Syntax:
%   r = raytrace(r0,surfaces)
%   r = raytrace(r0,surfaces,drawflag)
% r0 is a structure array of input rays, in the format defined in RAY.
% surfaces is a 1-by-ncomponents cell array, where each entry starts with a
%   function handle and continues with the necessary arguments to specify
%   the surface.  For example, for a planar surface followed by a
%   "spherical" (circular) surface in d=2, you'd have
%     surfaces{1} = {@opt2dline,lineparams}
%     surfaces{2} = {@opt2dcircle,circparams}
%   where the parameters are described in the appropriate file (opt2dline
%   or opt2dcircle).
% drawflag: if true, plots the rays as they are traced.
%
% Note that the surfaces have to be arranged in the order in which the
% ray will encounter them, otherwise weird things will happen...
% 
% See also: RAY.
  
  n_comp = length(comp);
  n_rays = length(r0);
  r = r0;
  if (nargin < 3)
    drawflag = false;
  end   
  
  for compIndex = 1:n_comp
    compdata = comp{compIndex};
    valid = [r.valid];
    validIndex = find(valid);
    [t,rnew] = feval(compdata{:},r(valid)); % pass only valid rays
    if drawflag
      % Plot on screen
      for rayIndex = 1:length(validIndex)
        x0 = r(validIndex(rayIndex)).x0;
        %x1 = rnew(rayIndex).x0;
        x1 = x0 + r(validIndex(rayIndex)).e * t(rayIndex);
        col = r(rayIndex).rgb;
        if (ndims(x0) == 3)
          line([x0(1) x1(1)],[x0(2) x1(2)],[x0(3) x1(3)],...
            'Color',col);
        else
          line([x0(1) x1(1)],[x0(2) x1(2)],...
            'Color',lightencolors(col,1-rnew(rayIndex).I));
        end
      end
    end
    r(valid) = rnew;
  end
 