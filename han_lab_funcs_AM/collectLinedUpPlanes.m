function collectLinedUpPlanes(nplanes)
    baseName1=input('saveName day1?');
    baseName2=input('saveName day2?');

for p=1:nplanes
    [fname1,pname1]=uigetfile('F*_proc.mat','Choose Suite2P Processed File For First Plane');
    [fname2,pname2]=uigetfile('F*_proc.mat','Choose Suite2P Processed File For Second Plane');
    [fnameReg,pnameReg]=uigetfile('*_reg.mat','Choose Suite2P Register File');

    day1=load([pname1,fname1]);
    day2=load([pname2,fname2]);
    load([pnameReg,fnameReg]); %Regi file. regi.dat.overlapping_rois has rois we want Nx2
    
    overRois=regi.dat.overlapping_rois;

    day1_inds=overRois(:,1);
    day2_inds=overRois(:,2);

    
    cellinds=find([day1.dat.stat.iscell]);
    F=day1.dat.Fcell{1}((day1_inds),:)';
    nF=day1.dat.FcellNeu{1}((day1_inds),:)';
    Fc=redo_dFF(F,day1.dat.ops.imageRate,30);
    Fc2=Fc;
    masks=zeros(length(day1_inds),size(day1.dat.res.iclust,1),size(day1.dat.res.iclust,2));
    for ii=1:length(day1_inds)
        masks(ii,:,:)=(day1.dat.res.iclust==day1_inds(ii));
    end    
    Fs=day1.dat.ops.imageRate;
    meanImage=squeeze(day1.dat.mimg(:,:,2));
    procImage=squeeze(day1.dat.mimg_proc(:,:,2));
    fName=day1.dat.filename;
    spks=day1.dat.sp{1}(day1_inds,:)';
    
    save([pname1,baseName1,'_SP_Registeredwith_',baseName2,num2str(p),'.mat'],'F','nF','Fc','Fc2',...
        'masks','Fs','meanImage','procImage','fName','spks')

    cellinds=find([day2.dat.stat.iscell]);
    F=day2.dat.Fcell{1}((day2_inds),:)';
    nF=day2.dat.FcellNeu{1}((day2_inds),:)';
    Fc=redo_dFF(F,day2.dat.ops.imageRate,30);
    Fc2=Fc;
    masks=zeros(length(day2_inds),size(day2.dat.res.iclust,1),size(day2.dat.res.iclust,2));
    for ii=1:length(day2_inds)
        masks(ii,:,:)=(day2.dat.res.iclust==day2_inds(ii));
    end    
    Fs=day2.dat.ops.imageRate;
    meanImage=squeeze(day2.dat.mimg(:,:,2));
    procImage=squeeze(day2.dat.mimg_proc(:,:,2));
    fName=day2.dat.filename;
    spks=day2.dat.sp{1}(day2_inds,:)';

    save([pname2,baseName2,'_SP_Registeredwith_',baseName1,num2str(p),'.mat'],'F','nF','Fc','Fc2',...
        'masks','Fs','meanImage','procImage','fName','spks')


end
    