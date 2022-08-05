function [fullcount, header] = olympus2imagine(basename, savedirpath)
%olympus2imagine takes Olympus files of the .oif format saved down by the ASV
%software on the FV1000 (support for other olympus models may be needed in future)
%and translates them into a header text file and a binary image data file that are
%functionally identical to those produced by the Imagine software (file formats and
%extensions will hereafter match those produced by Imagine). 'FV1000' will be
%indicated in the header.camera field to indicate the original data source. 
%
%It should be noted that you must first navigate to the outer directory where
%the properties '.oif' file is located for your experiment (data files are
%all saved inside a subdirectory but jumping in will be handled
%automatically by olympus2imagine)
% Syntax:
%   [fullcount, header] = olympus2imagine(basename);
%   [fullcount, header] = olympus2imagine(basename, savedirpath);
%   [fullcount, header] = olympus2imagine(basenames, savedirpath);
% where
%   "basename" is the original basename from the experimental data.
%      "basenames" is an optional syntax where a cell array of serial
%         names can be supplied for automatic concatenation (these
%         experiments need to have been collected under identical
%         circumstances). Present names in the order you want them
%         concatenated!  The final outcome will be named for the first
%         basename.
%   "savedirpath" is a string indicating the path to the directory to save
%              the translated files into.  If none is provided the file
%              will be saved into the present directory.
%   "fullcount" is the total number of pixel elements written to the binary
%              file and this number should be exactly equal to the product of image
%              height, width, and total frame count of the experiment.  
%   "header" is the name of the translated header
%  
%Copywrite 2006 Terrence Holekamp

homedir = cd;

if nargin < 2
    savedirpath = homedir;
end

if ~exist([savedirpath],'dir')
    [s] = mkdir(savedirpath);
    if ~s
        error('Cannot create output subdirectory. Check ''savedirpath'', write permissions, and disk free space')
    end
end

if ischar(basename)
    basename = cellstr(basename);
end
if ~iscell(basename)
    error('Experimental names are in an unrecogized or unsupported form')
end

savename = basename{1};
h = cell(1,length(basename));

for i = 1:length(basename)
    thisname = basename{i};
    cd(homedir);

    %load the FV1000 experiment properties text file as a header through
    %         imreadheader (which calls imreadheader_olympus)
    h{i} = imreadheader(thisname);

    %load serial tif images and save as binary file
    if strcmp(h{i}.scan_mode,'XYZT')
        if i == 1
            fullcount = 0;
        end
        for t = 1:h{i}.nstacks
            tindex = num2str(t);
            tlength = length(tindex);
            while tlength < 3
                tindex = ['0' tindex];
                tlength = length(tindex);
            end
            for z = 1:h{i}.frames_per_stack
                zindex = num2str(z);
                zlength = length(zindex);
                while zlength < 3
                    zindex = ['0' zindex];
                    zlength = length(zindex);
                end
                framename = ([thisname '_C001Z' zindex 'T' tindex '.tif']);
                cd([homedir filesep thisname '.oif.files']);
                if exist(framename,'file')
                    frame = imread(framename,'tif');
                else
                    frame = zeros(h{i}.width,h{i}.height,'uint16');
                end
                cd(savedirpath);
                if i == 1 && t == 1 && z == 1
                    fid = fopen([savename '.cam'], 'a');
                end
                count = fwrite(fid, frame, h{i}.prec);
                fullcount = fullcount + count;
            end
            fprintf('% d',t)
        end


    elseif strcmp(h{i}.scan_mode,'XYT')
        if i == 1
            fullcount = 0;
        end
        for t = 1:h{i}.nstacks
            tindex = num2str(t);
            tlength = length(tindex);
            while tlength < 3
                tindex = ['0' tindex];
                tlength = length(tindex);
            end
            framename = ([thisname '_C001T' tindex '.tif']);
            cd([homedir filesep thisname '.oif.files']);
            if exist(framename,'file')
                frame = imread(framename,'tif');
            else
                frame = zeros(h{i}.width,h{i}.height,'uint16');
            end
            cd(savedirpath);
            if i == 1 && t ==1
                fid = fopen([savename '.cam'], 'a');
            end
            count = fwrite(fid, frame, h{i}.prec);
            fullcount = fullcount + count;
            fprintf('% d',t)
        end

    else
        error('Scan mode unrecognized - unable to parse data file naming conventions')
    end
end

fclose(fid);

header = h{1};
oldframestring = ['FrameCount=' num2str(header.nframes)];
% oldstackstring = ['MaxSize=' num2str(header.nstacks)];%this line is probably not useful now that we use regexprep to fix the nstacks info in the original header
header.nframes = header.nframes*(length(basename));
header.nstacks = header.nstacks*(length(basename));
newframestring = ['FrameCount=' num2str(header.nframes)];
newstackstring = sprintf('$1MaxSize=%d\n$3', header.nstacks);

header.wholeheader = strrep(header.wholeheader,oldframestring,newframestring);
header.wholeheader = regexprep(header.wholeheader, '(\[Axis\ 4.*)(MaxSize=.*)(PixUnit=\"ms\")', newstackstring);

predictedfilesize = (header.height*header.width*header.nstacks*header.depth);

if ~isequal(fullcount,predictedfilesize)
    error('File write checksum failed - predicted file size did not match number of written elements');
end

save(savename, 'header');

cd(homedir);



