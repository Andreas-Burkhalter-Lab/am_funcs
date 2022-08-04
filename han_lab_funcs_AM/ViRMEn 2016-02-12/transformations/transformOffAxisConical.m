function coords2D = transformOffAxisConical(coords3D)

% define constants that describe setup geometry
S = 1280/768; % Aspect ratio of the projection (w/h)
T = 0.466476151587624; % Throw ratio of the projector (d/w)
H = 20; % Height of the screen
Rt = 11.5; % Top radius of the screen
Rb = 6; % Bottom radius of the screen
A = -2; % Location of the bottom ring with respect to the animal

% For 100% offset projector
gamma = atand(1/(2*S*T));

% Assume that the projector has been aligned to the top ring, make Xp a
% parameter and calculate Phi and Zp as a function of it
Xp = -4.5; %-10.5;
Phi = -atand(((Rt^2*S^2*T^2 + Rt^2 - Xp^2)^(1/2) - Rt*S*T)/(Rt - Xp));
Zp = H+A-(Rt+Xp)/tand(Phi);

% Offset direction (1 for x, 2 for y)
D = 2;

% Reused values
cs = cosd(Phi);
sn = sind(Phi);
B = A - H*Rb/(Rt-Rb);
M = H/(Rt-Rb);
Tv = 2*S*T;

% Parameters
a = B*Tv*cs;
b = -Tv*(Xp*cs + (Zp-B)*sn);
c = M*Tv*(Xp*cs+Zp*sn);
d = B*sn;
e = (Zp-B)*cs - Xp*sn;
f = -M*(Zp*cs-Xp*sn);
g = -Tv*tand(gamma);
h = B*Tv;

% create an output matrix of the same size as the input
coords2D = zeros(size(coords3D));

% by default, make all points visible
coords2D(3,:) = true;

% calculate radius
r = sqrt(coords3D(1,:).^2 + coords3D(2,:).^2);

% find points that do not project to the screen, lock them to the top or
% bottom ring, and set their visibility to 0
maxTan = (H+A)/Rt;
notOnScreen = coords3D(3,:)./r>maxTan;
coords3D(3,notOnScreen) = r(notOnScreen)*maxTan;
coords2D(3,notOnScreen) = false;
minTan = A/Rb;
notOnScreen = coords3D(3,:)./r<minTan;
coords3D(3,notOnScreen) = r(notOnScreen)*minTan;
coords2D(3,notOnScreen) = false;

% calculate transformations
coords2D(D,:) = (a*coords3D(D,:) + b*coords3D(3,:) + c*r) ./ (d*coords3D(D,:) + e*coords3D(3,:) + f*r) + g;
coords2D(3-D,:) = (h*coords3D(3-D,:)) ./ (d*coords3D(D,:) + e*coords3D(3,:) + f*r);