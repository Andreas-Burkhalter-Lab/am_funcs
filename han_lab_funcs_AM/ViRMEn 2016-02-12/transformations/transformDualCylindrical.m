function coords2D = transformDualCylindrical(coords3D)

s = 1;
p = 1;
scr = get(0,'screensize');

%original
aspectRatio = scr(3)/scr(4);
screenHeight = 21; %inches
radiusOfScreen = 11; %%inches
offsetX = 0;        %inches
offsetY = 0;  %inches
offsetZ = 0;

coords2D = zeros(size(coords3D));



%% Sid edit 7/26/17

coords2D(1,:) = coords3D(1,:)./coords3D(2,:);

disp(coords2D(1,1));   % diagnostic 

coords2D(1,:) = coords3D(1,:)./coords3D(2,:);
coords2D(2,:) = coords3D(3,:)./coords3D(2,:);

coords2D(3,:) = true;
theta = asin(coords2D(1,:)/aspectRatio);

realX = sin(theta)*radiusOfScreen + offsetX;
realY = (coords2D(2,:))*screenHeight/2 +offsetY*2/screenHeight;
realZ = cos(theta)*radiusOfScreen + offsetZ;

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
