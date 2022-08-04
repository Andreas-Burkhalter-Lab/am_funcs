% complete pipeline for calcium imaging data pre-processing
clear;
addpath(genpath('../NoRMCorre'));               % add the NoRMCorre motion correction package to MATLAB path
gcp;        % start a parallel engine
% foldername = 'G:\MA\AB Videos\161115_AB\';   
foldername = pwd;
        % folder where all the files are located. Currently supported .tif,
        % .hdf5, .raw, .avi, and .mat files
files = subdir(fullfile(foldername,'*plane*.mat'));   % list of filenames (will search all subdirectories)
FOV = [512,512];
numFiles = length(files);

%% motion correct (and save registered h5 files as 2d matrices (to be used in the end)..)
% register files one by one. use template obtained from file n to
% initialize template of file n + 1; 

non_rigid = true; % flag for non-rigid motion correction

template = [];
for i = 1:numFiles
    name = files(i).name;
    if non_rigid
        options_nonrigid = NoRMCorreSetParms('d1',512,'d2',512,'grid_size',[128,128],...
            'overlap_pre',64,'mot_uf',4,'bin_width',100,'max_shift',24,'max_dev',8,'us_fac',50,...
            'output_type','h5','h5_filename',[name(1:end-4),'_nr.h5']);
        [M,shifts,template] = normcorre_batch(name,options_nonrigid,template); 
        save([name(1:end-4),'_shifts_nr.mat'],'shifts','-v7.3');           % save shifts of each file at the respective subfolder
    else    % perform rigid motion correction (faster, could be less accurate)
        options_rigid = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'bin_width',100,'max_shift',32,...
            'output_type','h5','h5_filename',[name(1:end-4),'_rig.h5']);
        [M,shifts,template] = normcorre_batch(name,options_rigid,template); 
        save([name(1:end-4),'_shifts_rig.mat'],'shifts','-v7.3');           % save shifts of each file at the respective subfolder
    end
    
        video=uint16(h5read([name(1:end-4),'_nr.h5'],'/mov'));
        delete([name(1:end-4),'_nr.h5']); 
%         video=permute(squeeze(video), [2 1 3]);
        name=strrep(name,'plane','NC');
        save(name,'video','-v7.3');
    
end


