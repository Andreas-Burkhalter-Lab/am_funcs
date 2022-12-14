function [c,r] = spim_components(s)
% s is a structure with the following fields:
%   s.laserpos = amount to "de-focus" laser
%   s.cylR = radius of curvature of cylindrical lens
%   s.cylCT = center thickness of cylindrical lens
%   s.cylmat = string representing cylindrical lens material

  % Laser window
  xlw = -20;
  r = 3;
  p = struct('c',xlw-0.25,...
             'normal',[1 0],...
             'mat1','air',...
             'mat2','bk7',...
             'apc',[xlw-0.25 0],...
             'apr',r);
  c{1} = {@opt2dline,p};
  p = struct('c',xlw,...
             'normal',[1 0],...
             'mat1','bk7',...
             'mat2','air',...
             'apc',[xlw 0],...
             'apr',r);
  c{2} = {@opt2dline,p};
  

  % Thorlabs C220TM-A lens
  aparams = struct('curv',-.150884,'k',-.554366,'A',-4.3914e-6,...
                   'B',3.81456e-7,'C',0,'D',0);
  xca = -12;
  r = 3;
  p = struct('c',xca-5,...
             'normal',[1 0],...
             'mat1','air',...
             'mat2','c0550',...
             'apc',[xca-5 0],...
             'apr',r);
  c{3} = {@opt2dline,p};
  p = struct('xc',[xca 0],...
             'optax',[1 0],...
             'funcparams',aparams,...
             'xmax',r,...
             'mat1','c0550',...
             'mat2','air');
  c{4} = {@opt2daspheric,p};

  % The "cylindrical" surface of the cylindrical lens
  p = struct('xc',[0 0],...
             'Rvec',[-1 0]*s.cylR,...
             'theta',acos(1-s.cylCT/s.cylR),...
             'mat1',s.cylmat,...
             'mat2','air');
  c{5} = {@opt2dcircle,p};
  % The planar surface of the cylindrical lens
  xc = -s.cylR+s.cylCT;
  r = s.cylR*sin(p.theta);
  p = struct('c',xc,...
             'normal',[1 0],...
             'mat1',s.cylmat,...
             'mat2','saline',...
             'apc',[xc 0],...
             'apr',r);
  c{6} = {@opt2dline,p};
  % Now an extra planar surface to act as the "screen"
  screenpos = 9;
  p = struct('c',screenpos,...
             'normal',[1 0],...
             'mat1','saline',...
             'mat2','saline',...
             'apc',[screenpos 0],...
             'apr',3);
  c{7} = {@opt2dline,p};
  
  % Create the input rays
  r0.x0 = [xca-(5+7.96) + s.lasertwiddle,0];
  r0.w = 633;
  r0.e = [0 0];
  theta = linspace(-0.1,0.1,21);
  r = r0;
  for i = 1:length(theta)
    r(i) = r0;
    r(i).e = [cos(theta(i)) sin(theta(i))];
  end