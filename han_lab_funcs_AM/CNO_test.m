drive='F:/MA Data/Videos/';
% name=struct('folder',{'103015','110215','110415', '110915', '111315','112415','111615','112015'}, ...
%     'nums', {[3 4 5], [1 2 3], [1 2 4], [3 4 5], [0 1 2], [1 2 3],[1 2 3],[1 2 3]},...
%     'testname', {'M19 Saline 1', 'M19 CNO 1', 'M19 Saline 2', 'M20 Saline 1', 'M20 CNO 1', ...
%     'M20 CNO 2', 'M21 Saline 2', 'M21 CNO 1'});
% name=struct('folder',{'103015','110215','110415','111215', '110915', '111315','112415','111615','112015'}, ...
%     'nums', {[3 4 5], [1 2 3], [1 2 4], [3 4 5], [3 4 5], [0 1 2], [1 2 3],[1 2 3],[1 2 3]},...
%     'testname', {'M19 Saline 1', 'M19 CNO 1', 'M19 Saline 2', 'M19 CNO 2', 'M20 Saline 1', 'M20 CNO 1', ...
%     'M20 CNO 2', 'M21 Saline 2', 'M21 CNO 1'});
name=struct('folder',{'160307_MA','160307_MA','160308_MA','160308_MA'}, ...
    'nums',{[0 1 2],[3,5,6],[0 1 2],[3,5,6]},...
    'testname',{'M32 CNO1','M34 CNO1','M32 Saline 1','M34 Saline 1'});
meanF(4,3,4)=0;
% for testnum=1:2
%     filepath=[drive,name(testnum).folder,'/'];
%     for runnum=1:3
%        filename=[name(testnum).folder(1:end-2),'MA_000_00',num2str(name(testnum).nums(runnum))];
%        loadVideo(filepath,[filename,'.sbx'],4);
%     end
% end


    saveDir='F:/MA Data/CNO_F/';
    if ~exist([saveDir,],'dir')
        mkdir([saveDir]);
    end
for testnum=1:4
    filepath=[drive,name(testnum).folder,'/'];
    for runnum=1:3
        filename=[name(testnum).folder(1:end-2),'MA_000_00',num2str(name(testnum).nums(runnum))];
        for planenum=1:4
            load([filepath,filename,'_plane',num2str(planenum)]);
            imtool3D(chone(:,:,1:100));
            meanF(testnum,runnum,planenum)=mean(chone(:));
        end   
%         load([filepath,filename,'_plane',num2str(4),'_roibyhand__allF']);
%         meanF(testnum,runnum)=mean(chone(:));
    end
    label=name(testnum).testname;
    save([filepath,label,'_meanF.mat'],'meanF','label');
    meanFplanes = mean(meanF,3);
    figure;
    bar(meanFplanes(testnum,:));
%     figure
%     bar(meanF(testnum,:))
    title(label);
    saveas(gcf,[saveDir,label,'.jpg']);  
end

