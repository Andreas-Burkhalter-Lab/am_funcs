function [c,r] = objective_aspheric(s)
% s is a structure with the following fields:
%   s.wd = amount to "de-focus" fiber (translate the source)
%   s.immerse = string representing the material ('air' or 'saline') in
%     the back surface of the cylindrical lens is immersed.
%   s.lambda = wavelength (in nm)

  % Thorlabs C220TM-A lens
  %aparams = struct('curv',-.150884,'k',-.554366,'A',-4.3914e-6,...
  %                 'B',3.81456e-7,'C',0,'D',0);
  % Thorlabs 35033-A (C330TM?) lens: f=3.10, 0.68NA, AR coated for
  % 400-600nm. Lens has a "design wavelength" of 830nm, might be good for
  % IR microscopy. Diameter is 6.34mm, with a clear aperture of 4.3mm. In
  % air, used as a single infinity-corrected lens, the focus spot is
  % 2.1mm in front of the lens. In water, this might be 2.8mm (less, if
  % one is using a 2nd objective for the infinity focus). This yields an
  % absolute maximum working angle of 41 degrees.
  aparams1 = struct('curv',-.313623,'k',-12.66386,'A',-1.245834e-2,...
                   'B',3.711945e-3,'C',-5.122391e-4,'D',3.108578e-5);
  aparams2 = struct('curv',0.363636,'k',-.613916,'A',-5.88919e-4,...
                   'B',1.76602e-5,'C',-1.01025e-5,'D',3.91487e-6);
  % Thorlabs A414-A. This is a cylindrical-shape lens (might be easier to
  % mount), 4.5mm diameter. While it's designed to be inf-focused with a
  % BK7 laser window, let's say roughly 2mm to focus from front surface,
  % which with water correction yields something like an absolute maximum
  % working angle of 50 degrees. The lens looks like it might have a
  % weird shape, though.
  
  % Thorlabs 350350-A. Also cylindrical shape, 4.715mm diameter, efl
  % 4.5mm. Inf focus in air, 2.46mm in front of lens surface (note: I find
  % 2.25 at 510nm). Max working
  % angle 54 degrees, 3.3mm. This one seems to feature a wide range of
  % diffraction-limited performance (430-1550nm), suggesting it may be
  % tolerant of trouble. Exit aperature 3.7mm, air NA of 0.41. Implies
  % water NA 0.54 even before we work inside the focal length.
  % This looks like a good choice!
  % Note front surface is _almost_ plano (convex with small curvature)
  % Glass is CO550 (note may be written c0550 in opt_refrindex)
  % Another note: perhaps this lens, by itself, might make a good 30x
  % objective! efl 4.5mm*1.34 = efl 6mm, which implies 30x.
  aparams1 = struct('curv',0.052258,'k',0,'A',0,...
                   'B',0,'C',0,'D',0);
  aparams2 = struct('curv',-.347249,'k',-.6415948,'A',-3.148028e-4,...
                   'B',2.546471e-5,'C',2.814358e-6,'D',3.307336e-7);
  xmax = 2.05;
  glass = 'c0550';

  p = struct('xc',[0 0],...
             'optax',[1 0],...
             'funcparams',aparams1,...
             'xmax',xmax,...
             'mat1',s.immerse,...
             'mat2',glass);
  c{1} = {@opt2daspheric,p};
  p = struct('xc',[3.655 0],...
             'optax',[1 0],...
             'funcparams',aparams2,...
             'xmax',xmax,...
             'mat1',glass,...
             'mat2','air');
  c{2} = {@opt2daspheric,p};

  %ctmp = pcx([3 0],[-1 0],20,,s.cylmat,'air','air');
  
  % Now an extra planar surface to act as the "screen"
  screenpos = 50;
  p = struct('c',screenpos,...
             'normal',[1 0],...
             'mat1','air',...
             'mat2','air',...
             'apc',[screenpos 0],...
             'apr',screenpos/2);
  c{3} = {@opt2dline,p};
  
  % Create the input rays
  r0.x0 = [-s.wd,0];
  r0.w = s.lambda;
  r0.e = [0 0];
  r0.I = 1;
  sintheta = s.na/opt_refrindx(glass,s.lambda);
  thetamax = asin(sintheta);
  theta = linspace(-thetamax,thetamax,41);
  r = r0;
  for i = 1:length(theta)
    r(i) = r0;
    r(i).e = [cos(theta(i)) sin(theta(i))];
    r(i).I = exp(-theta(i)^2/(2*(s.na)^2));
  end
  