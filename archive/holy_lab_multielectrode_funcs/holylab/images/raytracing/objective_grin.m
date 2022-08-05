function [c,r] = objective_grin(s)
% [c,r] = objective_grin(s)
% where s is a structure with the following fields:
%   immerse: the immersion fluid at the GRIN entrance ('air' or 'saline')
%   wd: working distance
%   y: displacement from optic axis
%   pcx_pos: position of PCX lens

% GRIN lens
pitch = 0.18;
diam = 1.8;
grin_g = 0.332;
n0 = 1.5986;
R1 = inf;

% pitch = 0.11;
% diam = 3.0;
% grin_g = .203;
% n0 = 1.6249;
% R1 = inf;

%pitch = 0.2;
%diam = 1.8;
%grin_g = 0.423;
%n0 = 1.646;
%R1 = inf;

grinlength = 2*pi*pitch/grin_g;
grinparams = struct('xc',[0 0],'orient',[grinlength 0],'diam',diam,'n0', ...
		    n0,'g',grin_g,'R1',R1,'R2',inf,'mat1',s.immerse, ...
		    'mat2','air');

c{1} = {@opt2dgrin,grinparams};


% PCX lens
%ctmp = pcx([s.pcx_pos 0],[-1 0],30.9,4.9,'bk7','air','air');
%c{end+1} = ctmp{2};
%c{end+1} = ctmp{1};

% Set up a projection "screen"
screenpos = 20;
p = struct('c',screenpos,...
  'normal',[1 0],...
  'mat1','air',...
  'mat2','air',...
  'apc',[screenpos 0],...
  'apr',screenpos/2);
%  'apr',c{end}{2}.apr);
c{end+1} = {@opt2dline,p};



% Set up the rays
if ischar(s.wd)
  nm = opt_refrindx(s.immerse,s.lambda);
  C = (grinparams.n0 - nm)/(grinparams.n0*grinparams.R1);
  z = grin_g*grinlength;
  efl = nm / (grinparams.n0 * (grin_g*sin(z) + C*cos(z)));
  s.wd = efl * cos(z);
  disp(['Working distance: ' num2str(s.wd)])
  disp(['EFL: ' num2str(efl)]);
end
radius = grinparams.diam/2;
thetamin = atan((-radius-s.y)/s.wd);
thetamax = atan((radius-s.y)/s.wd);
theta = linspace(thetamin,thetamax,81);
r0.x0 = [-s.wd,s.y];
r0.e = [0 0];
r0.I = 1;
r0.w = s.lambda;
for i = 1:length(theta)
  r(i) = r0;
  r(i).e = [cos(theta(i)), sin(theta(i))];
end

