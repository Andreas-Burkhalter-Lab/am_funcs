%% TwoEnv abf
cellplanes(planes)=0;
for p=1:planes
    if manual==1;
        disp(['Plane ', num2str(p)]);
        [names{p},paths{p}]=uigetfile('*.mat','pick your files');
    else
        names{p}=varargin{4}{p};
        paths{p}=varargin{5}{p};
    end
    %     files{p}=name{p}(1:end-4);
    load([paths{p},names{p}]);
%     lp=designfilt('lowpassfir','FilterOrder',128, ...
%         'PassbandFrequency',5,'StopbandFrequency',7,...
%         'SampleRate',length(forwardvel)/length(F)*Fs);
    F0=F;
    F=Fc2;
    cellplanes(p)=size(F,2);
    %     F=bsxfun(@rdivide,F,max(F,[],1));
    %     F=bsxfun(@minus,F,mean(F,1));
    %Downsample ybinned to match F
    if length(ybinned) ~= size(F,1)
        rewratio=length(ybinned)/size(F,1);
        ratio=floor(length(ybinned)/size(F,1));
        path=downsample(ybinned(1:(ratio*size(F,1))),ratio);
        
        fornan=forwardvel;
        fornan(abs(fornan)<.02)=NaN;
        rotnan=rotationvel;
        rotnan(abs(rotnan)<.02)=NaN;
        tidx=1:length(forwardvel);
        forwardvel=interp1(tidx(~isnan(fornan)),fornan(~isnan(fornan)),tidx,'pchip','extrap');
        rotationvel=interp1(tidx(~isnan(rotnan)),rotnan(~isnan(rotnan)),tidx,'pchip','extrap');
        for jj=1:length(F)
            if (jj*rewratio)<length(ybinned)
                rotationvel(jj)=max(rotationvel(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
                forwardvel(jj)=max(forwardvel(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
                ybinned(jj)=max(ybinned(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
            else
                rotationvel(jj)=max(rotationvel(round(((jj-1)*rewratio+1)):end));
                forwardvel(jj)=max(forwardvel(round(((jj-1)*rewratio+1)):end));
                ybinned(jj)=max(ybinned(round(((jj-1)*rewratio+1)):end));
            end
        end
        rotationvel((jj+1):end)=[];
        forwardvel((jj+1):end)=[];
        forwardvel=smooth(forwardvel);
        rotationvel=smooth(rotationvel);
        ybinned((jj+1):end)=[];
    end
    
    %Split into novel and familiar
    
    novel_start=round(timeSplit(1)/rewratio);
    novel_end=round(timeSplit(2)/rewratio);
    
    envinds{1}=1:novel_start;
    envinds{2}=(novel_start+1):novel_end;
    envinds{3}=(novel_end+1):length(F);
    for env=1:length(env_label)
        forward{env}=forwardvel(envinds{env});
        rotation{env}=rotationvel(envinds{env});
        ybin{env}=ybinned(envinds{env});
        paths{env}=path(envinds{env})-min(path);
        paths{env}=paths{env}/max(paths{env})*180/5;
        pos{env}=ceil(paths{env}+eps);
        times{env}=linspace(1,length(ybin{env})/Fs,length(ybin{env}));
        
        if ~isempty(Fall{env})
            Fall{env}=[Fall{env},F(envinds{env},:)];
        else
            Fall{env}=F(envinds{env},:);
        end
        if ~isempty(F0all{env})
            F0all{env}=[F0all{env},F0(envinds{env},:)];
        else
            F0all{env}=F0(envinds{env},:);
        end
    end
end


%% SingleEnv abf

for p=1:planes
    if manual==1;
        disp(['Plane ', num2str(p)]);
        [names{p},paths{p}]=uigetfile('*.mat','pick your files');
    else
        names{p}=varargin{4}{p};
        paths{p}=varargin{5}{p};
    end
    %     files{p}=name{p}(1:end-4);
    load([paths{p},names{p}]);
%     lp=designfilt('lowpassfir','FilterOrder',128, ...
%         'PassbandFrequency',5,'StopbandFrequency',7,...
%         'SampleRate',10);
    F0=F;
    F=Fc2;
    %     Fsx=1:(1000/Fs):(length(F)*(1000/Fs));
    %     datax=1:length(ybinned);
    %     F=zscore(F);
    %     F=bsxfun(@rdivide,F,max(F,[],1));
    %     F=bsxfun(@minus,F,mean(F,1));
    %Downsample ybinned to match F
    if length(ybinned) ~= size(F,1)
        rewratio=length(ybinned)/size(F,1);
        fornan=forwardvel;
        fornan(abs(fornan)<.02)=NaN;
        rotnan=rotationvel;
        rotnan(abs(rotnan)<.02)=NaN;
        tidx=1:length(forwardvel);
        forwardvel=interp1(tidx(~isnan(fornan)),fornan(~isnan(fornan)),tidx,'pchip','extrap');
        rotationvel=interp1(tidx(~isnan(rotnan)),rotnan(~isnan(rotnan)),tidx,'pchip','extrap');
        for jj=1:length(F)
            if (jj*rewratio)<length(ybinned)
                rotationvel(jj)=max(rotationvel(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
                forwardvel(jj)=max(forwardvel(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
                ybinned(jj)=max(ybinned(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
            else
                rotationvel(jj)=max(rotationvel(round(((jj-1)*rewratio+1)):end));
                forwardvel(jj)=max(forwardvel(round(((jj-1)*rewratio+1)):end));
                ybinned(jj)=max(ybinned(round(((jj-1)*rewratio+1)):end));
            end
        end
        rotationvel((jj+1):end)=[];
        forwardvel((jj+1):end)=[];
        forwardvel=smooth(forwardvel);
        rotationvel=smooth(rotationvel);
        ybinned((jj+1):end)=[];
    end
%     forwardvel=smooth(forwardvel);
%     rotationvel=smooth(rotationvel);
    path=smooth(ybinned);
    path=(path-min(path));
    %     path=path/max(path)*1e6;
    pos=ceil(path+eps);
    if ~isempty(Fall)
        if length(Fall)<length(F)
            Fall=[Fall,F(1:length(Fall),:)];
        else
            Fall=[Fall(1:length(F),:), F];
        end
    else
        Fall=F;
    end
    %     novel_start=round(timeSplit(1)/rewratio);
    %     novel_end=round(timeSplit(2)/rewratio);
end

