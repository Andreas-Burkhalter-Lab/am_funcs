function convertSuite2Pdata_AM
% AM edited 2020/6/3.... adapted from convertSuite2Pdata by Moi Arriaga
%%%% only selection files from one directory at once
plotOpt = 0; % show plots
% alpha = scalefactor for dFF neuropil before subtracting from dFF roi...
%   ...0.7 alpha value taken from Kamani et al. 2016: 'Cooperative subnetworks of molecularly-similar interneurons in mouse neocortex' 
neuropil_alpha_val = 0.7; 
prefix_str = 'dff_';
keepspikes = 0; 

[fname,pname]=uigetfile('F*_proc.mat','Choose Suite2P Processed File','MultiSelect', 'on');
if ~iscell(fname)
    fname = {fname};
end
% if ~iscell(pname)
%     pname = {pname};
% end
nplanes = length(fname);

for p=1:nplanes
    load([pname,fname{p}])
    fprintf(['Getting df/f for ' fname{p} '...\n'])
    dat.stat = struct2table(dat.stat);
    dat.stat.iscell = logical(dat.stat.iscell);
    cellinds = find(dat.stat.iscell);
    F=dat.Fcell{1}((cellinds),:)';
    nF=dat.FcellNeu{1}((cellinds),:)';
    isnanF = find(any(isnan(F)));
    if ~isempty(isnanF) % check for empty traces, plot where they are in mean image
        isNanInd1 = cellinds(isnanF(1));
        figure
        imagesc(dat.mimg(:,:,2))
        colormap('gray')
        hold on
        scatter(dat.stat.xpix{isNanInd1}, dat.stat.ypix{isNanInd1},'r','.')
        error(['Found NaNs in F traces for ROIs [' num2str(isnanF) '].'])
    end
    
    % AM: Fc = get dFF by taking 8th percentile within (forward and backward) window as baseline
    dff=redo_dFF_AM(F,dat.ops.imageRate,30,nF,plotOpt,neuropil_alpha_val); 
%     Fc2=Fc; %%% AM commented out because Fc2 is same as Fc
    masks=zeros(length(cellinds),size(dat.res.iclust,1),size(dat.res.iclust,2));
    for ii=1:length(cellinds)
        masks(ii,:,:)=(dat.res.iclust==cellinds(ii));
    end    
    Fs=dat.ops.imageRate;
    meanImage=squeeze(dat.mimg(:,:,2));
    procImage=squeeze(dat.mimg_proc(:,:,2));
    fName=dat.filename;
    if isfield(dat,'sp') && keepspikes
        spks=dat.sp{1}(cellinds,:)';
    else 
        spks=[];
    end
    
    save([pname,prefix_str,strrep(fname{p},'F_','')],'F','nF','dff',... %'Fc2',... %%% AM commented out Fc2 because Fc2 is same as Fc
        'masks','Fs','meanImage','procImage','fName','spks','neuropil_alpha_val')

end
    

%% check for nans