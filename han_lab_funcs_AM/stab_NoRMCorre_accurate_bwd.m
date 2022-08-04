clear all
% function stab_NoRMCorre(non_rigid,
% complete pipeline for calcium imaging data pre-processing
% clear;

% addpath(genpath('../NoRMCorre'));               % add the NoRMCorre motion correction package to MATLAB path
gcp;        % start a parallel engine
% foldername = 'G:\MA\AB Videos\161115_AB\';
foldername = pwd;
% folder where all the files are located. Currently supported .tif,
% .hdf5, .raw, .avi, and .mat files
num_files=input('Enter number of files to run:');
% num_files=1;
files{num_files}=0;
filesb{num_files}=0;
paths{num_files}=0;
num_planes(num_files)=0;
for f=1:num_files
    [files{f},paths{f}]=uigetfile('.mat','pick your files');
%     files{f}=name(1:end-4);
end
for f=1:num_files
    load([paths{f},files{f}])
    if exist('video','var')
       chone=video;
       clear video;
    end
    chone=chone(:,:,length(chone):-1:1);
    filesb{f}=strrep(files{f},'.mat','b.mat');
    save([paths{f},filesb{f}],'chone','-v7.3');    
end
FOV = [512,601];
numFiles = length(files);

%% motion correct (and save registered h5 files as 2d matrices (to be used in the end)..)
% register files one by one. use template obtained from file n to
% initialize template of file n + 1;

non_rigid = true; % flag for non-rigid motion correction

template = [];
for i = 1:num_files
% for i=1
    name = filesb{i};
    if non_rigid
%         options_nonrigid = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'grid_size',[48,48],... %grid_size,[128,128]
%             'overlap_pre',32,'mot_uf',8,'bin_width',10,'max_shift',24,'max_dev',8,'us_fac',50,...
%             'output_type','h5','h5_filename',[name(1:end-4),'_nr.h5']);
        options_nonrigid = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'grid_size',[256,256],... %grid_size,[128,128]
            'overlap_pre',64,'mot_uf',6,'bin_width',100,'max_shift',24,'max_dev',8,'us_fac',50,...
            'output_type','h5','h5_filename',[name(1:end-4),'_nr.h5']);
        [M,shifts,template] = normcorre_batch(name,options_nonrigid,template);
        save(strrep(name,'plane','shifts_NC'),'shifts','-v7.3');           % save shifts of each file at the respective subfolder
        video=uint16(h5read([name(1:end-4),'_nr.h5'],'/mov'));
        delete([name(1:end-4),'_nr.h5']);        
    else    % perform rigid motion correction (faster, could be less accurate)
        options_rigid = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'bin_width',100,'max_shift',32,...
            'output_type','h5','h5_filename',[name(1:end-4),'_rig.h5']);
        [M,shifts,template] = normcorre_batch(name,options_rigid,template);
        save(strrep(name,'plane','rigShifts'),'shifts','-v7.3');           % save shifts of each file at the respective subfolder
        video=uint16(h5read([name(1:end-4),'_rig.h5'],'/mov'));
        delete([name(1:end-4),'_rig.h5']);
    end
    
    
    %         video=permute(squeeze(video), [2 1 3]);
    name=strrep(name,'plane','NCaccnew');
    save(name,'video','-v7.3');
    
end

files
files{1}
