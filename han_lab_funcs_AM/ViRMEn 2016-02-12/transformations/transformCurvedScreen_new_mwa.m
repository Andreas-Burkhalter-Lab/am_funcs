function coords2D = transformCurvedScreen_new_mwa(coords3D)
% %modified to work with two optoma 750ML projectors on screen with radius of
% %~12". keystone correction and 4:3 aspect ratio
% %screensize is 3840 x 1080 but get screensize only reads in one monitor
% %need to set pixels manually
% 
% 
% s = 1;
% %p = 1;
% % p = 1.1;
% p = 1;%1.1 seems to help with wobbly edge effects. cuts off edge
% 
% % scr = get(0,'screensize');
% scr_width = 3840; %total screen width in pix
% % scr_height = 1080; 
% scr_height = 1080; 
% scr = [1 1 scr_width scr_height];
% %originals
% % aspectRatio = scr(3)/scr(4);
% % screenLength = 55.292;  %%INCHES
% % screenHeight = 12.9066; %inches
% % distanceFromScreen = 14; %inches
% %  radiusOfScreen = 48; %%inches
% % offsetX = 0;        %inches
% % offsetY = 0;  %inches
% 
% 
% %for samsung curved
% aspectRatio = scr(3)/scr(4);
% %screenLength = 31.25;  %%INCHES true measure
% % screenLength = 50.25;  %%INCHES.  lokks better, increased aspect
% screenLength = 31.25;  %%INCHES.  lokks better, increased aspect
% 
% % screenHeight = 13.25; %inches
% % % screenHeight = 24; %inches.  gets wider with larger values until flips at 24" with distance 12".  
% % distanceFromScreen = 12; %inches
% % % radiusOfScreen = 48; %%inches
% % % radiusOfScreen = 118.11; %%inches
% % radiusOfScreen = 12; %%inches
% % %radiusOfScreen = 10; %%inches
% % offsetX = 0;        %inches
% % % offsetY = 0.5;  %inches
% % offsetY = 0;  %inches
% % 
% % 
% % 
% % coords2D = zeros(size(coords3D));
% % 
% % coords2D(1,:) = coords3D(1,:)./coords3D(2,:);
% % coords2D(2,:) = coords3D(3,:)./coords3D(2,:);
% % coords2D(3,:) = true;
% % realX = coords2D(1,:)*(screenLength/2)/aspectRatio + offsetX;
% % %realX = coords2D(1,:)*(screenLength)/aspectRatio + offsetX;%deleting the /2 helps with strange elongation effect when turning at wall
% % realY = (coords2D(2,:))*(screenHeight/2) +(offsetY*2)/screenHeight;
% % %realZ = distanceFromScreen -(radiusOfScreen.^2 - realX.^2).^.5 + radiusOfScreen;
% % realZ = distanceFromScreen -(radiusOfScreen.^2 - realX.^2).^.5 + radiusOfScreen;
% 
P=26; %inches
R=12.25; %inches
Py=21.5;
Px=20;
% 
% sy=(P-R)./(cos(asin(coords3D(3,:)/(P-R))));
% sx=sqrt((Px-abs(coords3D(1,:))).^2+(Py-coords3D(2,:)).^2);
% 
% % 
% % coords2D(1,:) = realZ.*coords2D(1,:)-realX;
% % coords2D(2,:) = realZ.*coords2D(2,:)-realY;
% 
% coords2D(1,:) = sx .* coords3D(1,:);
% 
% coords2D(2,:) = sy .* coords3D(2,:);
% 
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



% function coords2D = transformCylindrical(coords3D)

viewAngle = 90*pi/180; % how far beyond 90 degrees to the left and right the screen is behind the mouse
backwardAngle = viewAngle - pi/2;

% scr = getScreenSize();%%get(0,'screensize');
% scr = get(0,'screensize');
scr_width = 3840; %total screen width in pix
% scr_height = 1080; 
scr_height = 1080; 
scr = [1 1 scr_width scr_height];
aspectRatio = scr(3)/scr(4);
coords2D = zeros(size(coords3D));

% calculate horizontal angle theta for each point. theta begins at -0 and goes clockwise to +2viewAngle
theta = atan2(coords3D(1,:),coords3D(2,:))+viewAngle;
theta(theta<0)=nan;
theta(theta>(2*viewAngle))=nan;
theta=theta-pi/2;
% calculate screen x for each point
% coords2D(1,:) = aspectRatio./(viewAngle).*theta - aspectRatio;

% disp(['Coords 3D X min: ',num2str(min(coords3D(1,:))),' max: ', num2str(max(coords3D(1,:)))])
% disp(['Theta X min: ',num2str(min(theta)),' max: ', num2str(max(theta))])
% disp(['Aspect Ratio 3D X : ',num2str(aspectRatio),' View Angle: ', num2str(viewAngle)])

% min((((Px-R*sin(abs(theta-pi/2))).^2+(Py-R*cos(abs(theta-pi/2))).^2).^.5))
% max((((Px-R*sin(abs(theta-pi/2))).^2+(Py-R*cos(abs(theta-pi/2))).^2).^.5))
%Decent at .8 .4 .71
%Slight shallow bulgee 1 .7 .6    2.5 .7 .25
%.8 .4 .61

FudgeDisc=1; %Center Slope
Fudge2=.4; %Edge Slope
FudgeExp=1;
temp=((Px-FudgeDisc*R*sin(abs(theta))).^2+(Py-Fudge2*R*cos(abs(theta))).^2).^.5;
% temp=theta-pi/2;
% temp=((100-20*sin(abs(theta-pi/2))).^2+(100-20*cos(abs(theta-pi/2))).^2).^.5;
temp=temp-min(temp);
temp=-temp/max(temp)/2.+1; %2.8
% temp=temp.*((theta>pi/2)-(theta<pi/2));
% figure; plot(temp); title(['F1: ',num2str(FudgeDisc),' F2: ',num2str(Fudge2),' F3: ', num2str(Fudge3)])


coords2D(1,:) = theta.*temp*aspectRatio;
% coords2D(1,:) = theta*aspectRatio;

% disp([' X min: ',num2str(min((coords3D(1,:)))),' max: ', num2str(max((coords3D(1,:))))])

% disp(['Coords 2D X min: ',num2str(min(coords2D(1,:))),' max: ', num2str(max(coords2D(1,:)))])
% 
% calculate the distance in the xy plane to each point
r=(coords3D(1,:).^2 + coords3D(2,:).^2).^.5;
% r=((Px-FudgeDisc*R*sin(abs(theta))).^2+(Py-Fudge2*R*cos(abs(theta))).^2).^.5;
phi=atan2(coords3D(3,:),r);
% const=1;
% calculate the y screen coordinate for each point
coords2D(2,:) = coords3D(3,:)./r*3;
C2D=coords2D(2,:);
C2D=(C2D-min(C2D))/max((C2D-min(C2D)));
C2D=(C2D+1)/2;
% t2=temp;
% t2=((-temp+max(temp))/max(-temp+max(temp))+9)/10;
linearFromCenterX=((aspectRatio-abs(coords2D(1,:))))/(aspectRatio)*.2+.8; %0:1:0
sinFromCenter=(1-abs(sin(theta))).^2+(1.1-abs(cos(theta))).^2;
linearFromCenterY=((1-(coords2D(2,:)))/2)*.1+.9;
% squishFactor=(1+(linearFromCenterX.*linearFromCenterY)/2);
squishFactor= (sinFromCenter/max(sinFromCenter)*.25+.75).*linearFromCenterX;
coords2D(2,:) =  coords2D(2,:) ./ squishFactor;
% coords2D(2,:) = coords3D(3,:)./(r*.05).*((-1*temp+1.66).*(normpdf(coords3D(3,:),0,1)+.6));
% coords2D(2,find(abs(coords2D(2,:)>2)))=nan;
% disp(['Coords 2D initial Z min: ',num2str(min(coords2D(2,:))),' max: ', num2str(max(coords2D(2,:)))])

% temp=((Px-FudgeDisc*R*sin(abs(theta))).^2+(Py-Fudge2*R*cos(abs(theta))).^2).^.5;
% % temp=theta-pi/2;
% % temp=((100-20*sin(abs(theta-pi/2))).^2+(100-20*cos(abs(theta-pi/2))).^2).^.5;
% temp=temp-min(temp);
% temp=-temp/max(temp)/5.3+1;

coords2D(2,find(abs(coords2D(2,:))>2.5))=nan;
% disp(['Coords 2D denaned Z min: ',num2str(min(coords2D(2,:))),' max: ', num2str(max(coords2D(2,:)))])

mC2D=coords2D(2,:)-min(coords2D(2,:));
mC2D=mC2D/max(mC2D)/1.5;
mC2D=(1-mC2D/5);
 temp=temp/2+.5;
% coords2D(2,:)=1.*((coords2D(2,:)).*(temp).^FudgeExp);%.*cos(.25*(coords2D(2,:)+1));%.*(1-2*normpdf(coords2D(2,:)+1,0,1.5));

% coords2D(2,:)=coords2D(2,:)./(1-(coords2D(2,:)-min(coords2D(2,:)))/max(coords2D(2,:)-min(coords2D(2,:))));

% coords2D(2,:)=(coords2D(2,:)+1).^.9-1;
% coords2D(2,:)=5*coords2D(2,:).*(normpdf(coords2D(2,:)+1,1,3));

% disp(['Z scaling min: ',num2str(min(mC2D)),' max: ', num2str(max(mC2D))])
% disp(['Coords 3D Z min: ',num2str(min(coords3D(3,:))),' max: ', num2str(max(coords3D(3,:)))])
% disp(['r min: ',num2str(min(r)),' max: ', num2str(max(r))])
% disp(['Coords 2D Z min: ',num2str(min(coords2D(2,:))),' max: ', num2str(max(coords2D(2,:)))])
% 
% disp(' ')
% disp(' ')
% disp(' ')
% set all points to true
coords2D(3,:) = true;

% clip stuff left and right of the screen
f = find(abs(coords2D(1,:))>aspectRatio);
coords2D(1,f) = sign(coords2D(1,f))*aspectRatio;
coords2D(3,f) = false;

% clip stuff above and below the screen
f = find(abs(coords2D(2,:))>1);
coords2D(2,f) = sign(coords2D(2,f));
coords2D(3,f) = false;
