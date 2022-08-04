function coords2D = transformCurvedScreen(coords3D)
%foreshortening effect
%also elongation when turning at wall
s = 1;
%p = 1;
p = 1.1;
% p = 1.5;%1.1 seems to help with wobbly edge effects. cuts off edge

scr = get(0,'screensize');

%originals
% aspectRatio = scr(3)/scr(4);
% screenLength = 55.292;  %%INCHES
% screenHeight = 12.9066; %inches
% distanceFromScreen = 14; %inches
%  radiusOfScreen = 48; %%inches
% offsetX = 0;        %inches
% offsetY = 0;  %inches


%for samsung curved
aspectRatio = scr(3)/scr(4);
%screenLength = 31.25;  %%INCHES true measure
screenLength = 50.25;  %%INCHES.  lokks better, increased aspect
% screenLength = 55.25;  %%INCHES.  lokks better, increased aspect

screenHeight = 13.25; %inches
distanceFromScreen = 14; %inches
% radiusOfScreen = 48; %%inches
radiusOfScreen = 118.11; %%inches
% radiusOfScreen = 13; %%inches
% radiusOfScreen = 500.11; %%inches
%radiusOfScreen = 10; %%inches
offsetX = 0;        %inches
offsetY = 0.5;  %inches
% offsetY = 1;  %inches



coords2D = zeros(size(coords3D));

coords2D(1,:) = coords3D(1,:)./coords3D(2,:);
coords2D(2,:) = coords3D(3,:)./coords3D(2,:);
coords2D(3,:) = true;
realX = coords2D(1,:)*(screenLength/2)/aspectRatio + offsetX;
%realX = coords2D(1,:)*(screenLength)/aspectRatio + offsetX;%deleting the /2 helps with strange elongation effect when turning at wall
realY = (coords2D(2,:))*(screenHeight/2) +(offsetY*2)/screenHeight;
%realZ = distanceFromScreen -(radiusOfScreen.^2 - realX.^2).^.5 + radiusOfScreen;
realZ = distanceFromScreen -(radiusOfScreen.^2 - realX.^2).^.5 + radiusOfScreen;

coords2D(1,:) = realZ.*coords2D(1,:)-realX;
coords2D(2,:) = realZ.*coords2D(2,:)-realY;


f = coords3D(2,:)<=0;
coords2D(3,f) = false;
coords2D(1,f) = p*sign(coords3D(1,f))*aspectRatio;
coords2D(2,f) = p*sign(coords3D(3,f));

f = find(abs(coords2D(1,:))>aspectRatio);
coords2D(1,f) = p*sign(coords2D(1,f))*aspectRatio;
coords2D(3,f) = false;

f = find(abs(coords2D(2,:))>1);
coords2D(2,f) = p*sign(coords2D(2,f));
coords2D(3,f) = false;
