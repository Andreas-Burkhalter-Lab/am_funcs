function convertSuite2Pdata(nplanes)
% AM edited 18/5/17 to call redo_dFF_AM rather than redo_dFF
for p=1:nplanes
    [fname{p},pname{p}]=uigetfile('F*_proc.mat','Choose Suite2P Processed File');
end

baseName=input('saveName?');

for p=1:nplanes
        load([pname{p},fname{p}]);

    cellinds=find([dat.stat.iscell]);
    F=dat.Fcell{1}((cellinds),:)';
    nF=dat.FcellNeu{1}((cellinds),:)';
    % AM: Fc = get dFF by taking 8th percentile within (forward and backward) window as baseline
    Fc=redo_dFF_AM(F,dat.ops.imageRate,30,nF); % 1/22/18 AM added 4th arg nF..... Fc = 'F corrected' - see Arriaga and Han 2017
%     Fc2=Fc; %%% AM commented out because Fc2 is same as Fc
    masks=zeros(length(cellinds),size(dat.res.iclust,1),size(dat.res.iclust,2));
    for ii=1:length(cellinds)
        masks(ii,:,:)=(dat.res.iclust==cellinds(ii));
    end    
    Fs=dat.ops.imageRate;
    meanImage=squeeze(dat.mimg(:,:,2));
    procImage=squeeze(dat.mimg_proc(:,:,2));
    fName=dat.filename;
    if isfield(dat,'sp')
    spks=dat.sp{1}(cellinds,:)';
    else 
        spks=[];
    end
    save([pname{p},baseName,'_SP_',num2str(p),'.mat'],'F','nF','Fc',... %'Fc2',... %%% AM commented out Fc2 because Fc2 is same as Fc
        'masks','Fs','meanImage','procImage','fName','spks')

end
    