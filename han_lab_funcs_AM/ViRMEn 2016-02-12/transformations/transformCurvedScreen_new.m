function coords2D = transformCurvedScreen_new(coords3D)
%modified to work with two optoma 750ML projectors on screen with radius of
%~12". keystone correction and 4:3 aspect ratio
%screensize is 3840 x 1080 but get screensize only reads in one monitor
%need to set pixels manually


s = 1;
%p = 1;
% p = 1.1;
p = 1;%1.1 seems to help with wobbly edge effects. cuts off edge

% scr = get(0,'screensize');
scr_width = 3840; %total screen width in pix
% scr_height = 1080; 
scr_height = 1080; 
scr = [1 1 scr_width scr_height];
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
% screenLength = 50.25;  %%INCHES.  lokks better, increased aspect
screenLength = 31.25;  %%INCHES.  lokks better, increased aspect

screenHeight = 13.25; %inches
distanceFromScreen = 12; %inches
% radiusOfScreen = 48; %%inches
% radiusOfScreen = 118.11; %%inches
radiusOfScreen = 12; %%inches
%radiusOfScreen = 10; %%inches
offsetX = 0;        %inches
% offsetY = 0.5;  %inches
offsetY = 0;  %inches



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
