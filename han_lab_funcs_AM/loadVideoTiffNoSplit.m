function loadVideoTiffNoSplit(varargin)
if nargin<2
    [filename,filepath]=uigetfile('*.sbx','Choose SBX file');
    dir='uni';
elseif nargin<3
filepath=varargin{1};
filename=varargin{2};
dir='uni';
else
filepath=varargin{1};
filename=varargin{2};
dir=varargin{3};    
end
% win=60*Fs;

cd (filepath); %set path
stripped_filename=regexprep(filename,'.sbx','');
z = sbxread(stripped_filename,1,1);
global info;
chone = sbxread(stripped_filename,0,info.max_idx+1);
chone = squeeze(chone);
chone = double(chone);
framenum=size(chone,3);

javaaddpath 'C:\Program Files\MATLAB\R2018a\java\mij.jar'
javaaddpath 'C:\Program Files\MATLAB\R2018a\java\ij-1.51w.jar'
MIJ.start;    %calls Fiji

%%
numframes = size(chone,3);
lims(1)=min(chone(:))
lims(2)=max(chone(:))
lenVid=3000;
for ii=1:ceil(length(chone)/lenVid)
% ii=1;
    if ii>9
        currfile=strcat(stripped_filename,'_x',num2str(ii),'.mat');
    else
        currfile=strcat(stripped_filename,'_',num2str(ii),'.mat');
    end
    if strcmp(dir,'uni') || strcmp(dir,'uni new')
    chtemp=chone(:,:,((ii-1)*lenVid+1):min(ii*lenVid,length(chone)));
    elseif strcmp(dir,'bi')
        chtemp=chone(:,110:721,((ii-1)*lenVid+1):min(ii*lenVid,length(chone)));
    elseif strcmp(dir,'bi new')
            chtemp=chone(:,110:709,((ii-1)*lenVid+1):min(ii*lenVid,length(chone)));
    end
%     chtemp=chone;
    imageJ_savefilename=strrep([filepath,'\',currfile(1:end-4),'.tif'],'\','\\'); %ImageJ needs double slash
    imageJ_savefilename=['path=[' imageJ_savefilename ']'];
    MIJ.createImage('chone_image', gray2ind(mat2gray(chtemp,double(lims)),round(ceil(lims(2))/2)), true);
    MIJ.run('Save', imageJ_savefilename);   %saves with defined filename
    MIJ.run('Close All');
end

clear chone;

