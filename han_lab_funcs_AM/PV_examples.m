% plot_figs_for_Rotation_talk
%select how many days
saveDir='F:\MA Data\InterneuronsLowerCheck\';

nfiles=input('How Many Files');
name=input('ID');
Ffile{nfiles}=1;
ybin{nfiles}=1;
forward{nfiles}=1;
rotation{nfiles}=1;
cpp{nfiles}=1;
masksfiles{nfiles}=1;
Fraw{nfiles}=1;
%select m-files
for n=1:nfiles
   for p=1:4
        disp(['Plane ', num2str(p)]);
        [names{p},paths{p}]=uigetfile('*.mat','pick your files');
   end 
   [Ffile{n},ybin{n},forward{n},rotation{n},cpp{n},masksfiles{n},Fraw{n},~]=align_running(paths,names,4);
   load([paths{1},names{1}]);
    rewratio=length(rewards)/(length(Ffile{1}));
        novel_start=round(timeSplit(1)/rewratio);
    novel_end=round(timeSplit(2)/rewratio);
    rewinds{n}{1}=1:timeSplit(1);
    rewinds{n}{2}=(timeSplit(1)+1):timeSplit(2);
    rewinds{n}{3}=(timeSplit(2)+1):length(rewards);
    rewinds{n}{4}=1:length(rewards);
    envinds{n}{1}=1:novel_start;
    envinds{n}{2}=(novel_start+1):novel_end;
    envinds{n}{3}=(novel_end+1):length(Ffile);
    envinds{n}{4}=1:length(Ffile);
    
end
    plotcells=[3 10 13];
    
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    for n=1:nfiles
    for ii=1:length(plotcells)
    subplot(nfiles,3,(n-1)*3+ii)
    hold on
    cellF=Ffile{n}(:,plotcells(ii));
    plot(cellF)
        line([envinds{n}{1}(end) envinds{n}{1}(end)], [min(cellF) max(cellF)],'color','r')
        line([envinds{n}{2}(end) envinds{n}{2}(end)], [min(cellF) max(cellF)],'color','r')
    end
    end
    suptitle('PV Cell Remapping')
    saveas(gca,[saveDir,'Sample Cells ',name,'.jpg']);
    savefig([saveDir,'Sample Cells ',name,'.fig']);
    
    
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    for n=1:nfiles
    for ii=1:length(plotcells)
    subplot(nfiles,3,(n-1)*3+ii)
    hold on
    cellF=Ffile{n}(:,plotcells(ii));
        mcell(1)=mean(cellF(envinds{n}{1}));
        mcell(2)=mean(cellF(envinds{n}{2}));
        mcell(3)=mean(cellF((envinds{n}{2}+1):end));
        
        semcell(1)=1.96*std(cellF(envinds{n}{1}))/sqrt(length(envinds{n}{1}));
        semcell(2)=1.96*std(cellF(envinds{n}{1}+1:envinds{n}{2}))/sqrt(length(envinds{n}{1}+1:envinds{n}{2}));
        semcell(3)=1.96*std(cellF((envinds{n}{2}+1):end))/sqrt(length(cellF((envinds{n}{2}+1):end)));
        plot(mcell,'k:')
        errorbar(mcell,semcell,'r.')
        set(gca,'XTick',[1 2 3],'XTickLabel',{'Familiar','Novel','Familiar'})
        if max(mcell)<.3
            ylim([-.1 .3])
        else
            ylim([-.1 max(1.1*max(mcell),0)])
        end
    end
    end
    suptitle('PV Cell Remapping')
    saveas(gca,[saveDir,'Sample Cells Bar ',name,'.jpg']);
    savefig([saveDir,'Sample Cells Bar',name,'.fig']);
%plot

%plot