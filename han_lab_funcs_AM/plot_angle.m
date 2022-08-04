function plot_angle(F,up,down,varargin)
if nargin>3
pname=varargin{1};
name=varargin{2};
end
if nargin>5
part=varargin{3};
end
Fup=F(up,:);
Fdo=F(down,:);

mFup=nanmean(Fup);
mFdo=nanmean(Fdo);
nanloc=isnan(mFup)|isnan(mFdo);
mFup(nanloc)=[];
mFdo(nanloc)=[];

angle=dot(mFdo,mFup)/(norm(mFdo)*norm(mFup));

figure;
line([0 1],[0 0],'color','k');
hold on
line([0 cos(angle)], [0 sin(angle)],'color','k')
ylim([-.1 1.1]); xlim([-.1 1.1])
title(['Angle: ',num2str(angle*180/pi)])

if nargin>3
    saveas(gca,[pname,name,'\',name,' DifAngle.jpg']);
        saveas(gca,[pname,name,'\',name,' DifAngle.svg']);

end

if nargin>5
    import mlreportgen.dom.*;
    I=Image([pname,name,'\',name,' DifAngle.jpg']);
    append(part,TableEntry(I));
end

end