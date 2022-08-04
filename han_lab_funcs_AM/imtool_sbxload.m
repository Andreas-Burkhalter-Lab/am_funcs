% close all
%%%% from Ed Han
clear all

load_100=0;
[tifffilename,tiffpath]=uigetfile('*.sbx','pick your sbx file');
%tic
cd (tiffpath); %set path
stripped_tifffilename=regexprep(tifffilename,'.sbx','');
z = sbxread(stripped_tifffilename,1,1);
 global info;
if load_100==1 
    frameload=1000;
else
    frameload=info.max_idx+1; %to load all frames
end

chone = sbxread(stripped_tifffilename,0,frameload);
chone = squeeze(chone);
framenum=size(chone,3);
testmean=zeros(1,framenum);
    for i=1:framenum
        testmean(i)=mean(mean(chone(:,:,i)));
    end
figure,plot(testmean);

   tool=imtool3D(chone(:,:,:));
%  tool1=imtool3D(chone(:,:,1:4:end)); figure,plot(testmean(1:4:end));
%  tool2=imtool3D(chone(:,:,2:4:end));figure,plot(testmean(2:4:end));
%  tool3=imtool3D(chone(:,:,3:4:end));figure,plot(testmean(3:4:end));
%  tool4=imtool3D(chone(:,:,4:4:end));figure,plot(testmean(4:4:end));

%  tool=imtool3D(chone(:,50:650,:));
%tool=imtool3D(newmov(:,:,1:100));%does not cut off deadband

% 
% tool=imtool3D(newmov(:,50:650,1:4:end));
% 
% figure,plot(testmean(3:4:end));
% tool=imtool3D(newmov(:,50:650,3:4:end));
% 
% figure,plot(testmean(4:4:end));
% tool=imtool3D(newmov(:,50:650,4:4:end));

% figure,plot(testmean(1:16:end));
% tool=imtool3D(newmov(:,50:650,1:16:end));

% figure,plot(testmean(8:16:end));
% tool=imtool3D(newmov(:,50:650,8:16:end));
