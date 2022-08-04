%% TODO:
%Calculate Mutual Information between cell and position
%Sanity check with artificial data
%Figure out why the response changes
% function R = remapping (MouseID, Fc3, ybinned, Fs)

% function remapping(MouseID, Run, lFile, saveDir)
% close all; clear all;
%% Load and declare
saveDir='F:\MA Data\Figs\';
epo=input('Number of epochs');
planes=input('Number of planes');
MouseID = input('Mouse ID and Day');
files{epo,planes}=0;
paths{epo,planes}=0;
Fstruct{epo,planes}=0;
pos{epo,planes}=0;
Features{epo,planes}=0;
Probs{epo,planes}=0;
for e=1:epo
    disp(['Environment ', num2str(e)]);
    for p=1:planes
        disp(['Plane ', num2str(p)]);
        [name,paths{e,p}]=uigetfile('*.mat','pick your files');
        files{e,p}=name(1:end-4);
        load([paths{e,p},name]);
        
        F=Fc3;
        %Downsample ybinned to match Fc3
        if length(ybinned) ~= size(F,1)
            ratio=floor(length(ybinned)/size(F,1));
            ybinned=downsample(ybinned(1:(ratio*size(F,1))),ratio);
        end
        path=(ybinned-min(ybinned));
        path=path/max(path)*180/5;
        pos{e,p}=ceil(path+eps);
        
        Fstruct{e,p}=F;
    end
end

if ~exist([saveDir,MouseID,'\'],'dir')
    mkdir([saveDir,MouseID,'\']);
end

num_locs=max(pos{1,1});
% num_cells=size(Fstruct{1,1},2);


%% Spike distribution for each environment
%Feature count sums activity of neurons at each location (loc,neuron)
for e=1:epo
    for p=1:planes
        tvsdf_temp=Fstruct{e,p};
        pos_temp=pos{e,p};
        
        [feature_count, ~] = calculate_feature_count(tvsdf_temp,pos_temp);
        [feature_prob, ~] = calculate_feature_count(tvsdf_temp>0,pos_temp);
        
        %Prior prob shows the fraction of time mouse is in each loc
        prior_prob = histcounts(pos_temp,max(pos_temp));
        prior_probnorm = prior_prob/sum(prior_prob);
        
        feature_count = tuning_curve_processing(feature_count, prior_prob);
        feature_prob = tuning_curve_processing(feature_prob, prior_prob);
        
        Features{e,p} = feature_count;
        Probs{e,p} = feature_prob;
    end
end
%% Mutual Information
% % 
% %MI = sigma ( Pi (Ri/R) log2(Ri/R)) (Markus 1994)
% %Pi = probability for occupancy
% %Ri = Mean (firing rate) F for bin R = mean (firing rate) F for cell
% %Using feature_count, this is corrected for prior prob and gaussian
% %smoothed already. Try both ways?
% 
% bMI(1000,num_cells)=0;
% zMIs(num_cells,num_environs)=0;
% for j = 1:num_environs
%     RiR = bsxfun(@rdivide, features(j).Count, mean(features(j).Count));
%     %Add epsilon to prevent log(0) problems
%     RiR = bsxfun(@plus, RiR, eps);
%     MI = bsxfun(@times, prior_prob', (RiR.*log2(RiR)));
%     MI = sum(MI);
%     
%     %Bootstrap MI, try 1000 different sequences
%     for i=1:1000
%         fcountboot=calculate_feature_count(Fc3(randperm(length(Fc3)),:),pos);
%         fcountboot=tuning_curve_processing(fcountboot,prior_prob);
%         RiR = bsxfun(@rdivide, fcountboot, mean(fcountboot));
%         %Add epsilon to prevent log(0) problems
%         RiR = RiR + eps;
%         bMI(i,:) = sum(bsxfun(@times, prior_prob', (RiR.*log2(RiR))));
%     end
%     zM=zscore([MI; bMI]);
%     zMIs(:,i)=zM(1,:)';
% end
%% Find out at what location each cells fires each time
%This records the location at which each cell is firing
%multiply position vector by each activity vector - zero out the positions
%that were inactive for that cell
%cell_fired_rising uses rising edge activity
%shape is [path, cells]
% cell_fired_rising=bsxfun(@times,pos',Fc3br);
% cell_fired=bsxfun(@times,pos',(Fc3>0));

% 
% %posL, negL are the number of blocks return [0 1 1 1 0 0 2 2 0 3]
% %pos neg are the indices
% [posL,pos_run,negL,neg_run]=split_paths(ybinned);
% 
% % activity organized by location and run number
% tuning_changes=cell(num_cells,2,max([posL negL]));
% 
% for n = 1:num_cells
%     for i = 1:max(posL)
%         slope=pos_rect(posL==i)';
%         activity=Fc3(posL==i,n);
%         [sorted,sortind]=sort(slope);
%         tuning_changes{n,1,i}=[sorted activity(sortind)];
%     end
%     for i = 1:max(negL)
%         slope=pos_rect(negL==i)';
%         activity=Fc3(negL==i,n);
%         [sorted,sortind]=sort(slope);
%         tuning_changes{n,2,i}=[sorted activity(sortind)];
%     end
% end
% 
% % Plot tuning curves
% 
% if ~exist([saveDir,MouseID], 'dir')
%     mkdir([saveDir,MouseID]);
% end
% 
% %positive
% for n=1:num_cells
%     figure('units','normalized', 'Position', [.01 .05 .98 .87])
%     hold on
%     %     n=117;
%     for i = 1:max(posL)
%         temp=tuning_changes{n,1,i};
%         stem3(i*ones(size(temp(:,1))),temp(:,1),temp(:,2),'MarkerSize', 1)
%     end
%     title(['Positive running activity for neuron ',num2str(n), ' zMIstd=', num2str(zM(1,n))]);
%     view(-80,20)
%     h=gca;
%     saveas(h,[saveDir,'\',MouseID,'\', MouseID, Run, ' Cell ',num2str(n), ' Positive Runs.jpg']);
% end
% close all
% %negative
% for n=1:num_cells
%     figure('units','normalized', 'Position', [.01 .05 .98 .87])
%     hold on
%     for i = 1:max(negL)
%         temp=tuning_changes{n,2,i};
%         stem3(i*ones(size(temp(:,1))),temp(:,1),temp(:,2),'MarkerSize', 1)
%     end
%     title(['Negative running activity for neuron ',num2str(n), ' zMI=', num2str(zM(1,n))]);
%     view(-80,20)
%     h=gca;
%     saveas(h,[saveDir,'\',MouseID,'\', MouseID, Run, ' Cell ', num2str(n), ' Negative Runs.jpg']);
% 
% end
% close all
%% Color code on 2D stem

% colortime={'b','g','b'};
% tuning_changes = struct('env', []);
% for env=1:num_environs
%     [posL,pos_run,negL,neg_run]=split_paths(ybinned(indices_environ{env}));
%     pos_temp=pos_rect(indices_environ{env});
%     Fc3_temp=Fc3(indices_environ{env},:);
%     % activity organized by location and run number
%     tuning_changes(env).env=cell(num_cells,2,max([posL negL]));
%     for n = 1:num_cells
%         for i = 1:max(posL)
%             slope=pos_temp(posL==i)';
%             activity=Fc3_temp(posL==i,n);
%             [sorted,sortind]=sort(slope);
%             tuning_changes(env).env{n,1,i}=[sorted activity(sortind)];
%         end
%         for i = 1:max(negL)
%             slope=pos_temp(negL==i)';
%             activity=Fc3_temp(negL==i,n);
%             [sorted,sortind]=sort(slope);
%             tuning_changes(env).env{n,2,i}=[sorted activity(sortind)];
%         end
%     end
% end
% for n=1:num_cells
%     figure('units','normalized', 'Position', [.01 .05 .98 .87])
%     hold on
%     %     n=117;
%     for env=1:num_environs
%         for i = 1:size(tuning_changes(env).env,3)
%             temp=tuning_changes(env).env{n,1,i};
%             if ~isempty(temp)
%                 stem(temp(:,1),temp(:,2),'MarkerSize', 1,'color','b')
%             end
%         end
%     end
%     title(['Positive running activity for neuron ',num2str(n), ' zMIstd=', num2str(zM(1,n))]);
%     %         view(-80,20)
%     h=gca;
%     saveas(h,[saveDir,'\',MouseID,'\', 'Flat ', MouseID, Run, ' Cell ',num2str(n), ' Positive Runs.jpg']);
% end
% close all
% %negative
% for n=1:num_cells
%     figure('units','normalized', 'Position', [.01 .05 .98 .87])
%     hold on
%     for env=1:num_environs
%         for i = 1:size(tuning_changes(env).env,3)
%             temp=tuning_changes(env).env{n,2,i};
%             if ~isempty(temp)
%                 stem(temp(:,1),temp(:,2),'MarkerSize', 1,'color',colortime{env})
%             end
%         end
%     end
%     title(['Negative running activity for neuron ',num2str(n), ' zMI=', num2str(zM(1,n))]);
%     %         view(-80,20)
%     h=gca;
%     saveas(h,[saveDir,'\',MouseID,'\', 'Flat ', MouseID, Run, ' Cell ', num2str(n), ' Negative Runs.jpg']);
%     
% end
% close all


%}
%% Place cells
%
% %Preallocate
% place_cells_height=zeros(num_cells,1);
% place_cells_proportions=zeros(num_cells,1);
% numheight=0;
% numprop=0;
% winsize=3;
% for i=1:num_cells
% location_activity=feature_count(:,i);
% [maxval, maxloc] = max(location_activity);
% if maxval > (mean(location_activity) + std(location_activity))
%     numheight = numheight + 1;
%     place_cells_height(numheight)=i;
% end
% if maxloc - winsize <1
%     start = 1;
% else start = maxloc-winsize;
% end
% if maxloc + winsize > num_locs
%     stop = num_locs;
% else stop = maxloc + winsize;
% end
% if sum(location_activity((start):(stop))) / sum(location_activity)>.5
%    numprop = numprop + 1;
%    place_cells_proportions(numprop) = i;
%
% end
%
% end
% %Cleave off extra spaces
% place_cells_height = place_cells_height(1:numheight);
% place_cells_proportions = place_cells_proportions(1:numprop);
% %

% use=ones(num_cells,1);
%% Use cells with good MI
% tvsdfm=tvsdf(:,zM(1,:)>4);
% feature_countm=feature_count1(:,zM(1,:)>4);
% z=zM(zM(1,:)>4);
% zM=z;
% tvsdf=tvsdfm; feature_count1=feature_countm;
%% Plotting results


for p=1:planes
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    for e=1:epo
        subplot(1,epo,e)
        features=Features{e,p};
        [~,maxind] = max(features);
        [maxes, inds] = sort(maxind);
        imagesc(features(:,inds))
        colormap(jet)
        title(['Environment ', num2str(e),' plane ', num2str(p)]);
    end
end
% for i = 1:num_environs
%     subplot(2,num_environs,i);
%     imagesc(features(i).Count);
%     colormap(jet(4096));
%     colorbar;
%     title([name_environs{i}]);
%     subplot(2,num_environs,i+num_environs);
%     imagesc(features(i).Prob);
%     colorbar;
%     title([name_environs{i},' probability of firing'])
% end
% h=gcf;
% saveas(h,[saveDir,MouseID,'\', MouseID, Run, ' heatmaps.jpg']);
% 
% figure('units','normalized', 'Position', [.01 .05 .98 .87]);
% for env = 1:num_environs
%     temp=Fc3(indices_environ{env},:);
%     plot(temp(:,1))
%     subplot(2,num_environs,env);
%     nz = temp>0;
%     imagesc(-nz')
%     colormap('gray')
%     title(['Raster-like plot of activity in ', name_environs{env}])
%     freezeColors;
%     subplot(2,num_environs,env+num_environs);
%     plot(pos(indices_environ{env}))
%     
% end
% h=gcf;
% saveas(h,[saveDir,'\',MouseID,'\', MouseID, Run, ' rasters.jpg']);
% 
% % figure;
% % bar(zMIs)
% 
% R=0;
% % save(lFile,'features','-append');
% 
% figure('units','normalized', 'Position', [.01 .05 .98 .87]);
% title('Mean Flourescence by location');
% colormap(jet(4096*2));
% % use_neurons=[14 3 29 32 10 7 4 9 16 36 21 2 1 33 11 22 8 16]; mouse 1
% use_neurons=[4 31 52 40 13 8 1 10 11 9 15 22 50 36 44 2 27 49 ];% mouse 3
% % use_neurons=[22 16 6 8 25 4 27 14 11 9 10]; mouse 3.2
% % use_neurons=[107 8 46 93 131 103 44];
% heatmap=features(1).Count(:,use_neurons); 
% smheatmap=imresize(heatmap, [length(heatmap)*1000, length(use_neurons)]);
% for n=1:length(use_neurons)
%    subplot(1,length(use_neurons),n)
%    imagesc(smheatmap(:,n));
%    set(gca, 'YTickLabel',[],'XTickLabel',[]);
% %    title(['Sample ', num2str(n)]);
%   
% end
% set(gcf, 'Title','Mean Flourescence by location');
% %  colorbar; 
% time=linspace(1,length(Fc3)/Fs,length(Fc3));
% selected=1100:2250;
% % selected=1:length(tvsdf);
% tselected=1:length(selected);
% figure('units','normalized', 'Position', [.01 .05 .98 .87]);
% subplot(2,10,[11:18])
% hold on
% plot(time(tselected),pos(selected)/max(pos),'color','k')
% % line([0 10], [0 0],'color','k','LineWidth',2);
% xlabel('Seconds');
% ylabel('Position');
% % set(gca, 'YTickLabel',[]);
% 
% colors='br';
% for n=1:2
%     n=1;
%     hold on
%     subplot(2,10,[1:8])
%     plot(time(tselected),Fc2(selected,use_neurons(n))/max(Fc2(selected,use_neurons(n)))+1.05,colors(n))
%     if n==1 
%         xloc=1;
%     else xloc=length(tselected)/Fs;
%     end
% %     line([xloc xloc], [1.05 (1.05+(round(max(Fc2(selected,use_neurons(n))))/max(Fc2(selected,use_neurons(n))))/2)],'color',colors(n),'LineWidth',2)
%     max(Fc2(selected,use_neurons(n)))
%     round(max(Fc2(selected,use_neurons(n))))
%     round(max(Fc2(selected,use_neurons(n))))/max(Fc2(selected,use_neurons(n)));
%     xlabel('Seconds')
%     ylabel('dF/F')
%     set(gca, );
% %     axis off
% end
% subplot(2,10,[19])
% imagesc(smheatmap(:,1));
% colormap(jet(4096*2));
% set(gca, 'YTickLabel',[],'XTickLabel',[]);
% ylabel('Position')
% 
% subplot(2,10,[20])
% imagesc(smheatmap(:,2));
% colormap(jet(4096*2));
% set(gca, 'YTickLabel',[],'XTickLabel',[]);
% 
% % set(gca, 'YTickLabel',[],'XTickLabel',[]);
% % ylabel('position')
% 
% % 
% % end