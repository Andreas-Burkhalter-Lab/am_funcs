%runVideos
%Allows you to select multiple sbx files, separates out and stabilizes each
%plane, currently using the HMM method from SIMA
clear all; close all;
num_files=input('Enter number of files to run:');
% num_files=1;
files{num_files}=0;
paths{num_files}=0;
num_planes(num_files)=0;
for f=1:num_files
    [name,paths{f}]=uigetfile('*.sbx','pick your files');
    files{f}=name(1:end-4);
    num_planes(f)=input('Enter number of planes for this file: ');
end

for i=1:num_files
    %Run python stabilization code
    cd(paths{i});
    file=['"',paths{i},files{i},'"'];
    %To use hmm
    systemCommand = ['python E:\MA\Python_Scripts\simahmm.py ', file, ' ', num2str(num_planes(i))];
    system(systemCommand)
    %Reads in hdf5 file format from python code and converts it back into
    %mat files
    for p=1:num_planes(i)
        video=uint16(h5read([files{i},'_hmm',num2str(p),'.hdf5'],'/imaging'));
        delete([files{i},'_hmm',num2str(p),'.hdf5']); 

        video=permute(squeeze(video), [2 1 3]);
        save([paths{i},files{i},'_hmm',num2str(p),'.mat'],'video','-v7.3');
    end
    
 end
