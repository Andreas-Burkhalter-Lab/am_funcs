function loadVideoTiffNoSplit_AM(varargin)

% AM edited 19/7/22 on thermaltake


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
if nargin>3
    is_sbx = varargin{4};
end
if nargin>4
    file_ext = varargin{5};
end
% win=60*Fs;

cd (filepath); %set path
stripped_filename=regexprep(filename,'.sbx','');
if is_sbx
    z = sbxread(stripped_filename,1,1);
end
global info;


%% AM edited section for large files..... 4/8/18 on thermaltake
%    computer can't handle chunks larger than 30gb
if fsize([filename file_ext]) > 3e4 % if larger than 30GB
    if ~is_sbx
        error('edit this section so it works on non-sbx')
    end
    z = squeeze(z);
    framenum = info.max_idx+1;
    chone = NaN(size(z,1),size(z,2),framenum); % initialize at full size
    nFramesAtOnce = 1000;
    nchunks = ceil([info.max_idx+2]/nFramesAtOnce);
    for ichunk = 1:nchunks 
        ichunk
        startframe = (ichunk-1)*nFramesAtOnce + 1; %%% counting by chone frames, not sbxread frames, so first frame==1, not 0
        nFramesThisChunk = min([nFramesAtOnce, framenum-[(ichunk-1)*nFramesAtOnce]]);
        endframe = startframe + nFramesThisChunk - 1; %%% counting by chone frames, not sbxread frames, so first frame==1, not 0
        chone(:,:,startframe:endframe) = double(squeeze(sbxread(stripped_filename,startframe-1,nFramesThisChunk))); %%% load in only a small number of frames at once
    end
else
    if is_sbx
        chone = sbxread(stripped_filename,0,info.max_idx+1);
        chone = squeeze(chone);
        chone = double(chone);
    else % if is tif
        full_fname = [filepath filename file_ext];
        fileinfo = imfinfo(full_fname);
        framenum = length(fileinfo);
        xsize = fileinfo(1).Width;
        ysize = fileinfo(1).Height;
        chone = NaN(ysize,xsize,framenum); % initialize
        for iframe = 1:framenum
            chone(:,:,iframe) = imread(full_fname,iframe);
        end
    end
end
%%

framenum=size(chone,3);

javaaddpath 'C:\Program Files\MATLAB\R2018a\java\mij.jar'
javaaddpath 'C:\Program Files\MATLAB\R2018a\java\ij-1.51w.jar'
MIJ.start;    %calls Fiji

%%
numframes = size(chone,3);
lims(1)=min(chone(:))
lims(2)=max(chone(:))
lenVid=3000; %% number of frames per tiff chunk
for ii=1:ceil(length(chone)/lenVid) %% ii is the index of tiff chunk to be made
% ii=1;
    if ii>9
        currfile=strcat(stripped_filename,'_x',num2str(ii),'.mat');
    else
        currfile=strcat(stripped_filename,'_',num2str(ii),'.mat');
    end
    if ~is_sbx || strcmp(dir,'uni') || strcmp(dir,'uni new') %%% if tif rather than sbx, don't crop the image
        chtemp=chone(:,:,((ii-1)*lenVid+1):min(ii*lenVid,size(chone,3)));
    elseif strcmp(dir,'bi')
        chtemp=chone(:,110:721,((ii-1)*lenVid+1):min(ii*lenVid,length(chone))); % maybe change length(chone) to size(chone,3)
    elseif strcmp(dir,'bi new')
        chtemp=chone(:,110:709,((ii-1)*lenVid+1):min(ii*lenVid,length(chone))); % maybe change length(chone) to size(chone,3)
    end
%     chtemp=chone;
    imageJ_savefilename=strrep([filepath,'\',currfile(1:end-4),'.tif'],'\','\\'); %ImageJ needs double slash
    imageJ_savefilename=['path=[' imageJ_savefilename ']'];
    MIJ.createImage('chone_image', gray2ind(mat2gray(chtemp,double(lims)),round(ceil(lims(2))/2)), true);
    MIJ.run('Save', imageJ_savefilename);   %saves with defined filename
    MIJ.run('Close All');
end

clear chone;

