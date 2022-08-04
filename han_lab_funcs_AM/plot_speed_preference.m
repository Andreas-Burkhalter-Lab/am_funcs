function plot_speed_preference(Sp,forward,rotation,varargin)
if nargin>3
pname=varargin{1};
name=varargin{2};
type=varargin{3};
else type=[];
end
if nargin>6
part=varargin{4};
end
%% Calc Preferences
corrRunning=nan(size(Sp,2),1);
total=sqrt(forward.^2+rotation.^2);
for ii=1:length(corrRunning)
   tc=corrcoef(total(:,ii),Sp(:,ii));
   corrRunning(ii)=tc(2,1);
end
[~,inds]=sort(corrRunning);


%% Speed Preference
figure('units','normalized', 'Position', [.01 .05 .95 .86]);
s1=subplot(2,10,1:10);
imagesc(Sp(:,inds)')
title(['Cell Activity Raster with ',type])
colormap(flipud(gray))
ylabel('Cells')
% for ii=1:length(startPup)
% patch([startPup(ii) startPup(ii) stopPup(ii) stopPup(ii)],[0 size(Sp,2) size(Sp,2) 0],...
% 'b','FaceAlpha',.1,'EdgeAlpha',0)
% end
% for ii=1:length(startPdo)
% patch([startPdo(ii) startPdo(ii) stopPdo(ii) stopPdo(ii)],[0 size(Sp,2) size(Sp,2) 0],...
% 'r','FaceAlpha',.1,'EdgeAlpha',0)
% end

s2=subplot(2,10,11:20);
plot(sqrt(forward(:,1).^2+rotation(:,1).^2),'k')
hold on
title('Ball Speed')
plot(forward(:,1),'b--','LineWidth',.2)
plot(rotation(:,1),'r--','LineWidth',.2)
legend({'Total','Forward','Rotation'})
xlabel('Frames');
ylabel('Running Speed')
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

saveas(gca,[pname,name,'\',name,' SpeedRaster',type,'.jpg']);
saveas(gca,[pname,name,'\',name,' SpeedRaster',type,'.svg']);


if nargin>5
    import mlreportgen.dom.*;
    I=Image([pname,name,'\',name,' SpeedRaster',type,'.jpg']);
    append(part,TableEntry(I));  
end