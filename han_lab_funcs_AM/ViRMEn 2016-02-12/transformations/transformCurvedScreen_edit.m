function coords2D = transformCurvedScreen_edit(coords3D)
%not actually transformCurvedScreen_edit
%actually "transformCylindrical" from Brad
%Could not get to work, gives error on line 323 of virmenEngine
%Probably need to list name of transform somewhere in code, unknown where
%pasted this code into existing functional transform
%original transformCurvedScreen_edit at bottom

%viewAngle = 112.5*pi/180; % how far beyond 90 degrees to the left and right the screen is behind the mouse
%viewAngle = 70*pi/180;
viewAngle = 100*pi/180;%samsung
backwardAngle = viewAngle - pi/2;
 

scr = getScreenSize();%%get(0,'screensize');
aspectRatio = scr(3)/scr(4);
aspectRatio2 = scr(3)/scr(4);%mwa

coords2D = zeros(size(coords3D));

% calculate horizontal angle theta for each point. theta begins at -0 and goes clockwise to +2viewAngle
theta = atan2(coords3D(1,:),coords3D(2,:))+viewAngle;

% calculate screen x for each point
coords2D(1,:) = aspectRatio2./(viewAngle).*theta - aspectRatio2;

% calculate the distance in the xy plane to each point
%r=(coords3D(1,:).^2 + coords3D(2,:).^2).^.5;
r=(coords3D(1,:).^2 + coords3D(2,:).^2).^.5;
% calculate the y screen coordinate for each point
coords2D(2,:) = coords3D(3,:)./r;

% set all points to true
coords2D(3,:) = true;

% clip stuff left and right of the screen
f = find(abs(coords2D(1,:))>aspectRatio);
coords2D(1,f) = sign(coords2D(1,f))*aspectRatio;

%below to cut off most peripheral screens
% f = find(abs(coords2D(1,:))>aspectRatio*.59);
% coords2D(1,f) = sign(coords2D(1,f))*aspectRatio*.59;

%below to keep only central screen
% f = find(abs(coords2D(1,:))>aspectRatio*.20);
% coords2D(1,f) = sign(coords2D(1,f))*aspectRatio*.20;

coords2D(3,f) = false;

% clip stuff above and below the screen
f = find(abs(coords2D(2,:))>1);
coords2D(2,f) = sign(coords2D(2,f));
coords2D(3,f) = false;

%below is original transformCurvedScreen_edit
% s = 1;
% p = 1;
% 
% %change made
% scr = getScreenSize();%%get(0,'screensize');
% %change made
% 
% aspectRatio = scr(3)/scr(4);
% %screenLength = 55.292;  %%INCHES
% screenLength = 59.5; %inches BAR
% %screenHeight = 12.9066; %inches
% screenHeight = 21; %inches
% %distanceFromScreen = 14; %inches
% %radiusOfScreen = 48; %%inches
% distanceFromScreen = 21.25; %inches
% radiusOfScreen = distanceFromScreen; %%inches
% 
% offsetX = 0;  %inches
% offsetY = 0;  %inches
% 
% coords2D = zeros(size(coords3D));
% 
% coords2D(1,:) = coords3D(1,:)./coords3D(2,:); % unitless Xworld/Zworld
% coords2D(2,:) = coords3D(3,:)./coords3D(2,:); % unitless Yworld/Zworld
% coords2D(3,:) = true;
% 
% % realX = coords2D(1,:) + offsetX.*aspectRatio/screenLength;
% 
% realX = coords2D(1,:).*screenLength/2/aspectRatio + offsetX; % inches
% realY = coords2D(2,:).*screenHeight/2 + offsetY; % inches
% %realZ = distanceFromScreen -(radiusOfScreen.^2 - realX.^2).^.5+ radiusOfScreen;
% realZ = radiusOfScreen./cos(asin(realX./radiusOfScreen)); % inches
% % 
% coords2D(1,:) = realZ.*coords2D(1,:) - realX;
% coords2D(2,:) = realZ.*coords2D(2,:) - realY;
% 
% % Xreal = Xworld*screenLength/2; % inches
% % Yreal = Yworld*screenHeight/2; % inches
% % Zreal = radiusOfScreen./cos(asin(Xreal./radiusOfScreen)); % inches
% 
% %        (inches/unitless*unitless - inches)*(1/inches)         
% % Xscreen = (Zreal./Zworld.*Xworld - Xreal)./screenLength;
% % Yscreen = (Zreal./Zworld.*Yworld - Yreal)./screenHeight;
% % 
% % coords2D(1,:) = Xscreen;
% % coords2D(2,:) = Yscreen;
% % coords2d(3,:) = true;
% 
% f = coords3D(2,:)<=0;
% coords2D(3,f) = false;
% coords2D(1,f) = p*sign(coords3D(1,f))*aspectRatio;
% coords2D(2,f) = p*sign(coords3D(3,f));
% 
% f = find(abs(coords2D(1,:))>aspectRatio);
% coords2D(1,f) = p*sign(coords2D(1,f))*aspectRatio;
% coords2D(3,f) = false;
% 
% f = find(abs(coords2D(2,:))>1);
% coords2D(2,f) = p*sign(coords2D(2,f));
% coords2D(3,f) = false;
