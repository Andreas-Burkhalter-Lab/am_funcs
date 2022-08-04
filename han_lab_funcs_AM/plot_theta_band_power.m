function plot_theta_band_power(varargin)
set(0,'DefaultFigureColormap',jet)

saveDir='G:\MA\MA Data\Interneurons\';
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
for p=1:planes
    if manual==1;
        disp(['Plane ', num2str(p)]);
        [names{p},paths{p}]=uigetfile('*.mat','pick your files');
    else
        names{p}=varargin{4}{p};
        paths{p}=varargin{5}{p};
    end
    %     files{p}=name{p}(1:end-4);
end
for p=1:planes
        load([paths{p},names{p}]);

if exist('chone_corr','var')
    video=chone_corr;
    clear chone_corr
end


video=uint16(bsxfun(@minus,single(video),single(mean(mean(video,1),2))));

rmsfull=rms(video,3).^2;
figure
imagesc(rmsfull)
if ~exist([saveDir,'\',MouseID,'\'],'dir')
    mkdir([saveDir,'\',MouseID,'\'])
end
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Full RMS squares', '.jpg']);
bp=zeros(size(video,1),size(video,2));
for i=1:size(video,1)
    bp(i,:)=bandpower(squeeze(single(video(i,:,:)))',Fs,[4 12]);
end
figure
imagesc(bp)
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Theta Band power','.jpg']);

end
end