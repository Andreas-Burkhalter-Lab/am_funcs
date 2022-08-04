function [Fall, rewards,ybinned]=getData(varargin)
saveDir='F:\MA Data\Interneurons\';
if nargin<1
    Fs=input('Fs? ');
    planes=input('Number of planes: ');
    MouseID = input('Mouse ID and Day: ');
    names{planes}=0;
    paths{planes}=0;
    manual=1;
else
    Fs=varargin{1};
    planes=varargin{2};
    MouseID=varargin{3};
    manual=0;
end
Fall=[];


for p=1:planes
    if manual==1;
        disp(['Plane ', num2str(p)]);
        [names{p},paths{p}]=uigetfile('*.mat','pick your files');
    else
        names{p}=varargin{4}{p};
        paths{p}=varargin{5}{p};
    end
    load([paths{p},names{p}]);
    lp=designfilt('lowpassfir','FilterOrder',128, ...
        'PassbandFrequency',5,'StopbandFrequency',7,...
        'SampleRate',length(forwardvel)/length(F)*Fs);
    F=Fc2;
    if length(ybinned) ~= size(F,1)
        rewratio=length(ybinned)/size(F,1);

        for jj=1:200:(length(forwardvel)-200)
            forwardvel(jj:(jj+199))=max(forwardvel(jj:(jj+199)));
            rotationvel(jj:(jj+199))=max(rotationvel(jj:(jj+199)));
        end
        forwardvel=filtfilt(lp,forwardvel);
        rotationvel=filtfilt(lp,rotationvel);
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
        ybinned((jj+1):end)=[];
    end
    path=smooth(ybinned);
    path=(path-min(path));
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
end



end