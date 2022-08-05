pcxlength = 4.9;
R = 30.9;
w = 488;
s.y = 0;  % offset from optic axis
s.wd = 10;
mat1 = 'saline';
mat2 = 'air';
matl = 'bk7';
cpcx = pcx([0 0],[-1 0],30.9,4.9,matl,mat2,mat1);
cpcx = cpcx([2 1]);

% Create a "GRIN" that is really a PCX lens (g = 0)
diam = 2*cpcx{2}{2}.apr;
n0 = opt_refrindx(matl,w);
grinparams = struct('xc',[-pcxlength 0],'orient',[pcxlength 0],'diam',diam,'n0', ...
		    n0,'g',1e-8,'R1',R,'R2',inf,'mat1',mat1, ...
		    'mat2',mat2);
cgrin = {@opt2dgrin,grinparams};

% Draw both lenses (one on top of the other)
cla
feval(cpcx{1}{:});
feval(cpcx{2}{:});
feval(cgrin{:});

% Now an extra planar surface to act as the "screen"
screenpos = 5;
p = struct('c',screenpos,...
  'normal',[1 0],...
  'mat1','air',...
  'mat2','air',...
  'apc',[screenpos 0],...
  'apr',30);
cscreen = {@opt2dline,p};


% Set up some rays
radius = diam/2;
thetamin = atan((-radius-s.y)/(s.wd-pcxlength));
thetamax = atan((radius-s.y)/(s.wd-pcxlength));
theta = linspace(thetamin,thetamax,81);
r0.x0 = [-s.wd,s.y];
r0.w = 488;
r0.I = 1;
r0.e = [0 0];
for i = 1:length(theta)
  r(i) = r0;
  r(i).e = [cos(theta(i)), sin(theta(i))];
end


% Trace these rays
for i = 1:length(r)
  raytrace(r(i),[cpcx,{cscreen}],[1 0 0]);
  %raytrace(r(i),[{cgrin},{cscreen}],[0 1 0]);
end

shg
