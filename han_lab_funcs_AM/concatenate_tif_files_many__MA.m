clear all;
close all;

%get all files names

% window=100;   %original default
window=100;
% cated_filename1='140306_004-005.tif'; num_files1=2;
%   cated_filename2='140306_008-009.tif'; num_files2=2;
%     cated_filename3='121211_001-004_3.tif'; num_files3=4;
%    cated_filename4='131202_005-008_1.tif'; num_files4=4;
%    cated_filename5='131125_006-009_2.tif'; num_files5=4;
%    cated_filename6='131125_006-009_3.tif'; num_files6=4;
%    cated_filename7='16_003-005_1.tif'; num_files7=3;
%    cated_filename8='16_003-005_2.tif'; num_files8=3;
%    cated_filename9='16_003-005_3.tif'; num_files9=3;
%  cated_filename10='130307_014-017_2.tif'; num_files10=4;
%  cated_filename11='130307_014-017_3.tif'; num_files11=4;

%Don't specify different files for planes now. Just the one name and plane
%number is appended
cated_filenames={'1118_MA_000_000-003','1119_MA_002-005'};
num_groups=length(cated_filenames);
% cated_filenames={'1119_MA_000_002-005'};
num_files=[4,4];
num_planes=[3,3];


%   cated_tiff_filename12='130306_020-023_2.tif'; num_files12=4;
% cated_tiff_filename13='130306_020-023_3'; num_files13=4;
% cated_tiff_filename14='121024_006-007_2'; num_files14=2;
% cated_tiff_filename15='121024_006-007_3'; num_files15=2;
% cated_tiff_filename16='121809_47_7to9'; num_files16=3;
% cated_tiff_filename17='121809_47_10to12'; num_files17=3;
% cated_tiff_filename18='122109_49_1to3'; num_files18=3;
% cated_tiff_filename19='122109_49_4to6'; num_files19=3;
% cated_tiff_filename20='122109_49_7to9'; num_files20=3;
% cated_tiff_filename21='122109_50_1to3'; num_files21=3;
% cated_tiff_filename22='122109_50_4to7'; num_files22=4;
% cated_tiff_filename23='122109_50_8to10'; num_files23=3;
% cated_tiff_filename24='122109_47_1to3'; num_files24=3;
% cated_tiff_filename25='122109_47_4to6'; num_files25=3;
% cated_tiff_filename26='122109_47_7to9'; num_files26=3;
% cated_tiff_filename27='122109_47_10to12'; num_files27=3;
% cated_tiff_filename28='122309_49_1to3'; num_files28=3;
% cated_tiff_filename29='122309_49_4to6'; num_files29=3;
% cated_tiff_filename30='122309_49_7to9'; num_files30=3;
% cated_tiff_filename31='122309_49_10to14'; num_files31=5;
% cated_tiff_filename32='122309_50_1to3'; num_files32=3;
% cated_tiff_filename33='122309_50_4to7'; num_files33=4;
% cated_tiff_filename34='122309_50_8to10'; num_files34=3;
% cated_tiff_filename35='122309_50_11to14'; num_files35=4;
% cated_tiff_filename36='122309_47_1to3'; num_files36=3;
% cated_tiff_filename37='122309_47_4to6'; num_files37=3;
% cated_tiff_filename38='122309_47_7to9'; num_files38=3;

%120815 EH
h=warndlg('Close open directory folders!');%need to close folders to write to files
waitfor(h);
pause(1);


saveloc=struct('path', [], 'file', struct('name',[]));
for k=1:num_groups
    %     eval(['num_files=num_files' num2str(k)]);
    for j=1:num_files(k)
        %         eval(['cated_tiff_filename' num2str(k)]);
        cated_filenames{k}
        [filename,fpath]=uigetfile('*.sbx','pick your sbx file');
        
        %         eval(['fullfilename' num2str(k) '_' num2str(j) ' = [fpath filename]'])
        saveloc(k).path=fpath;
        saveloc(k).file(j).name=filename;
        cd (fpath); %set path
        %                 %for moving all scripts to "Local"  120628  EH
    end
end

%cd (tiffpath); %set path
%for moving all scripts to "Local"  120628  EH
%%

for k=1:num_groups
    %     eval(['num_files=num_files' num2str(k)]);
    %     eval(['cated_filename=cated_filename' num2str(k)]);
    
    %open files and concatenate; scale for same instensities
    for p=1:num_planes(k)
        cated_movie=[];
        
        for j=1:num_files(k)
            
            %         eval(['fullfilename=fullfilename' num2str(k) '_' num2str(j)]);
            fullfilename=[saveloc(k).path saveloc(k).file(j).name]
            
            stripped_filename=regexprep(fullfilename,'.sbx','');
            z = sbxread(stripped_filename,1,1);
            global info;
            newmov = sbxread(stripped_filename,0,info.max_idx+1);
            % newmov = sbxread(stripped_tifffilename,0,1000);
            newmov = squeeze(newmov);
            if j==1
                newmov=newmov(:,:,1:6660);
            end
            framenum=size(newmov,3);
            M=size(newmov,1);
            N=size(newmov,2);
            %Original code doesnt seem to parse out different frames??
            %         info=imfinfo(fullfilename);
            %         numframes=length(info);
            %         M=info(1).Width;
            %         N=info(1).Height;
            %
            %         chone=zeros(N,M,numframes);
            %         for i=1:numframes
            %             if mod(i,100)==1
            %                 j
            %                 i
            %             end
            %             chone(:,:,i)=imread(fullfilename,i,'Info',info);
            %         end
            
            %scale movie for seamless intensities
            chone=newmov(:,:,p:num_planes(k):end);
            clear newmov;
            if j==1
                meanlastframes=median(mean(mean(chone(:,:,1:150))));
            end
            
            meanfirstframes=median(mean(mean(chone(:,:,1:150))));
            chone=chone*(meanlastframes/meanfirstframes);
            size(chone)
            size(cated_movie)
            cated_movie=cat(3,cated_movie,chone);
            size(cated_movie)
            meanlastframes=median(mean(mean(chone(:,:,end-150:end))));
        end
        junk=squeeze(mean(mean(cated_movie)));
        numframes=length(junk);
        mean_all=mean(mean(mean(cated_movie)));
        
        junk2=zeros(size(junk));
        for kk=1:length(junk)
            cut=junk(max(1,kk-window):min(numframes,kk+window));
            cutsort=sort(cut);
            a=round(length(cut)*.08);
            junk2(kk)=cutsort(a);
        end
        
        %     cated_movie_1=zeros(size(cated_movie));
        for i=1:numframes
            if mod(i,100)==1
%                 i
            end
            %         cated_movie_1(:,:,i)=(cated_movie(:,:,i)/junk2(i))*mean_all;
            cated_movie(:,:,i)=(cated_movie(:,:,i)/junk2(i))*mean_all;
            
        end
        %     cated_movie=cated_movie_1;
        %     clear cated_movie_1;
        %
        
        save([saveloc(k).path, cated_filenames{k}, '_plane', num2str(p)],'cated_movie','-v7.3');
        cated_movie=[];
        clear 'cated_movie';
    end
    
    %baseline subtract whole movie
    
    
    %   pause(5);
    %iowritemovie_tif(cated_movie,tiffpath,cated_tiff_filename);
    %   saveastiff(uint16(cated_movie),cated_tiff_filename);
    %clear cated_movie;
    
    %For saving tif with ImageJ (Fiji).
    %Update Fiji and add Fiji scripts folder to matlab path.
    %Java heap size in Matlab is limit, 8Gb for EBH comp.
    % %     final_filename=[fpath cated_filenames{k}]; %files need to have path
    % %     imageJ_savefilename=strrep(final_filename,'\','\\'); %ImageJ needs double slash
    % %     imageJ_savefilename=['path=[' imageJ_savefilename ']'];
    % %     Miji;    %calls Fiji
    % %     %uint8 works but think it clips high rather than scaling. noisy.
    % %     MIJ.createImage('chone_image', uint16(cated_movie), true); %creates ImageJ file with 'name', matlab variable name
    % %     MIJ.run('Save', imageJ_savefilename);   %saves with defined filename
    % %     MIJ.run('Close All');
    % %     MIJ.exit;
    % %     %   cd (fpath); %set path. Miji seems to reset path to default
    % %
    % %
    % %
    % %
    % %
    % %
    % %
    % %
    % %     pause(5);
    % %
end

% beep; pause(1);
% beep; pause(1);
% beep; pause(1);
% % g = fittype('a*exp(x)+b');
% [gg,GOODNESS,OUTPUT]=FIT(x',a,'exp2');
% plot(a);hold on; plot(gg)