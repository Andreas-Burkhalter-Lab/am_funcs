%% Plotting results


figure('units','normalized', 'Position', [.01 .05 .98 .87]);

for i = 1:num_environs
    subplot(2,num_environs,i);
    imagesc(features(i).Count);
    colormap(jet(4096));
    colorbar;
    title([name_environs{i}]);
    subplot(2,num_environs,i+num_environs);
    imagesc(features(i).Prob);
    colorbar;
    title([name_environs{i},' probability of firing'])
end
h=gcf;
saveas(h,[saveDir,MouseID,'\', MouseID, Run, ' heatmaps.jpg']);

figure('units','normalized', 'Position', [.01 .05 .98 .87]);
for env = 1:num_environs
    temp=tvsdf(indices_environ{env},:);
    plot(temp(:,1))
    subplot(2,num_environs,env);
    nz = temp>0;
    imagesc(-nz')
    colormap('gray')
    title(['Raster-like plot of activity in ', name_environs{env}])
    freezeColors;
    subplot(2,num_environs,env+num_environs);
    plot(pos(indices_environ{env}))
    
end
h=gcf;
saveas(h,[saveDir,'\',MouseID,'\', MouseID, Run, ' rasters.jpg']);

% figure;
% bar(zMIs)

R=0;
% save(lFile,'features','-append');

figure('units','normalized', 'Position', [.01 .05 .98 .87]);
title('Mean Flourescence by location');
colormap(jet(4096*2));
% use_neurons=[14 3 29 32 10 7 4 9 16 36 21 2 1 33 11 22 8 16]; mouse 1
use_neurons=[4 31 52 40 13 8 1 10 11 9 15 22 50 36 44 2 27 49 ];% mouse 3
% use_neurons=[22 16 6 8 25 4 27 14 11 9 10]; mouse 3.2
% use_neurons=[107 8 46 93 131 103 44];
heatmap=features(1).Count(:,use_neurons); 
smheatmap=imresize(heatmap, [length(heatmap)*1000, length(use_neurons)]);
for n=1:length(use_neurons)
   subplot(1,length(use_neurons),n)
   imagesc(smheatmap(:,n));
   set(gca, 'YTickLabel',[],'XTickLabel',[]);
%    title(['Sample ', num2str(n)]);
  
end
set(gcf, 'Title','Mean Flourescence by location');
%  colorbar; 
time=linspace(1,length(tvsdf)/Fs,length(tvsdf));
selected=1100:2250;
% selected=1:length(tvsdf);
tselected=1:length(selected);
figure('units','normalized', 'Position', [.01 .05 .98 .87]);
subplot(2,10,[11:18])
hold on
plot(time(tselected),pos(selected)/max(pos),'color','k')
% line([0 10], [0 0],'color','k','LineWidth',2);
xlabel('Seconds');
ylabel('Position');
% set(gca, 'YTickLabel',[]);

colors='br';
for n=1:2
    n=1;
    hold on
    subplot(2,10,[1:8])
    plot(time(tselected),Fc2(selected,use_neurons(n))/max(Fc2(selected,use_neurons(n)))+1.05,colors(n))
    if n==1 
        xloc=1;
    else xloc=length(tselected)/Fs;
    end
%     line([xloc xloc], [1.05 (1.05+(round(max(Fc2(selected,use_neurons(n))))/max(Fc2(selected,use_neurons(n))))/2)],'color',colors(n),'LineWidth',2)
    max(Fc2(selected,use_neurons(n)))
    round(max(Fc2(selected,use_neurons(n))))
    round(max(Fc2(selected,use_neurons(n))))/max(Fc2(selected,use_neurons(n)));
    xlabel('Seconds')
    ylabel('dF/F')
    set(gca, );
%     axis off
end
subplot(2,10,[19])
imagesc(smheatmap(:,1));
colormap(jet(4096*2));
set(gca, 'YTickLabel',[],'XTickLabel',[]);
ylabel('Position')

subplot(2,10,[20])
imagesc(smheatmap(:,2));
colormap(jet(4096*2));
set(gca, 'YTickLabel',[],'XTickLabel',[]);

% set(gca, 'YTickLabel',[],'XTickLabel',[]);
% ylabel('position')

% 
% end