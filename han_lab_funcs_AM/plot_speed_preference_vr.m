function plot_speed_preference_vr(Sp,forward,varargin)
if nargin>2
pname=varargin{1};
name=varargin{2};
type=varargin{3};
end
if nargin>5
part=varargin{4};
import mlreportgen.dom.*;
end
%% Calc Preferences
corrRunning=nan(size(Sp,2),1);
total=sqrt(forward.^2);
for ii=1:length(corrRunning)
   tc=corrcoef(total(:,1),Sp(:,ii));
   corrRunning(ii)=tc(2,1);
end
[~,inds]=sort(corrRunning);


%% Speed Preference
figure('units','normalized', 'Position', [.01 .05 .95 .86]);
ha=tight_subplot(2,1,[.01 .05],[.05 .01],[.03 .01]);
% axes(ha(1));
s1=subplot(2,10,1:10);
imagesc(Sp(:,inds)')
colormap(flipud(gray))
ylabel('Cell')
title('Cell Activity Raster Organized by VR Speed')
% for ii=1:length(startPup)
% patch([startPup(ii) startPup(ii) stopPup(ii) stopPup(ii)],[0 size(Sp,2) size(Sp,2) 0],...
% 'b','FaceAlpha',.1,'EdgeAlpha',0)
% end
% for ii=1:length(startPdo)
% patch([startPdo(ii) startPdo(ii) stopPdo(ii) stopPdo(ii)],[0 size(Sp,2) size(Sp,2) 0],...
% 'r','FaceAlpha',.1,'EdgeAlpha',0)
% end

s2=subplot(2,10,11:20);
% axes(ha(2));
plot(sqrt(forward(:,1).^2),'k')
hold on
title('VR Speed')
plot(forward(:,1),'b--','LineWidth',.2)
ylabel('Speed')
xlabel('Frames')
legend({'Abs','Directional'})
axis(gca,'tight')
% for ii=1:length(startPup)
% patch([startPup(ii) startPup(ii) stopPup(ii) stopPup(ii)],[0 3 3 0],...
% 'b','FaceAlpha',.1,'EdgeAlpha',0)
% end
% for ii=1:length(startPdo)
% patch([startPdo(ii) startPdo(ii) stopPdo(ii) stopPdo(ii)],[0 3 3 0],...
% 'r','FaceAlpha',.1,'EdgeAlpha',0)
% end
linkaxes([s1 s2],'x')
% linkaxes([ha],'x')
if nargin>2
saveas(gca,[pname,name,'\',name,' SpeedRaster',type,'.jpg']);
saveas(gca,[pname,name,'\',name,' SpeedRaster',type,'.svg']);

end

if nargin>5
    import mlreportgen.dom.*;
    I=Image([pname,name,'\',name,' SpeedRaster',type,'.jpg']);
    append(part,TableEntry(I));  
end