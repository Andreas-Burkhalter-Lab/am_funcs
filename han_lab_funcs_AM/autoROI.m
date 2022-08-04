% function autoROI(pathplane,fileplane,pathref,planeref)
% num_mice=input('Enter number of mice to ROI:');
% num_planes{num_mice}=0;
% for n=1:num_mice
%    num_days{n}=input(['Number of days for mouse ',num2str(n)]);
%    num_planes{n}=input(['Number of planes for mouse ',num2str(n)]);
%
% end
% num_files=1;
num_days=input(['Number of days for mouse ']);
num_planes=input(['Number of planes for mouse']);
fileplanes{num_days}=0;
pathplanes{num_planes}=0;
% num_planes(num_days)=0;
fileref{num_planes}=0;
pathref{num_planes}=0;
for p=1:num_planes
    [fileref{p},pathref{p}]=uigetfile(['Select Reference Frame for plane ',num2str(p)]);
    for d=1:num_days
        [fileplanes{(p-1)*num_planes+d},pathplanes{(p-1)*num_planes+d}]=uigetfile('*.mat',['pick your motion correct video for plane ',num2str(p),' day ',num2str(d)]);
    end
end

% [fileplane,pathplane]=uigetfile('Motion corrected plane to add rois');
% [fileref,pathref]=uigetfile('Select Reference Frame');
numwindow=250;
for p=1:num_planes
    maskref=load([pathref{p},fileref{p}],'masks');
    ref_centroids=find_centroids(maskref.masks);
    numref=size(maskref.masks,1);
    
    for d=1:num_days
        load([pathplanes{(p-1)*num_planes+d},fileplanes{(p-1)*num_planes+d}]);
        
        newmasks=zeros(size(maskref.masks,1),size(video,1),size(video,2));
        frame=mean(video,3);
        rois=bwlabel(frame>(mean(frame(:)+3*std(frame(:)))),4);
        rois2=bwareaopen(rois,125);
        rinfo=regionprops(rois2);
        num_rois=length(rinfo);
        centroids=cat(1,rinfo.Centroid);
        area=cat(1,rinfo.Area);
        distmat=pdist2(ref_centroids,centroids);
        topleft=find(distmat(1,:)==min(distmat(1,:)));
        nx=centroids((topleft),1);
        ny=centroids((topleft),2);
        distx1=nx-ref_centroids(1,1);
        disty1=ny-ref_centroids(1,2);
        for i=1:numref
            tempmask=zeros(size(frame));
%             m=squeeze(maskref.masks(i,:,:)); using handrawn ROI as
%             template
            m=rois(rois==(i));
            [r,c]=find(m>0);
            mean(r)
            mean(c)
            r=round(r+disty1);
            c=round(c+distx1);
            mean(r)
            mean(c)
            locs=sub2ind(size(frame),r,c);
            tempmask(locs)=ones(size(tempmask(locs)));
            if i>1
                m=tempmask;
                newinfo=regionprops(bwlabel(m));
                newcentroids=cat(1,newinfo.Centroid);
                newdists=pdist2(ref_centroids,newcentroids);
                closest=find(newdists==min(newdists));
                tempmask=zeros(size(frame));
                nx=newcentroids(1);
                ny=newcentroids(2);
                distx=nx-ref_centroids(closest,1);
                disty=ny-ref_centroids(closest,2);
                [r,c]=find(m>0);
                mean(r)
                mean(c)
                r=round(r+disty);
                c=round(c+distx);
                mean(r)
                mean(c)
                locs=sub2ind(size(frame),r,c);
                tempmask(locs)=ones(size(tempmask(locs)));
            end
            
            newmasks(i,:,:)=tempmask;
        end
        
        %%
        numframes=length(video);
        N=size(video,2);
        M=size(video,3);
        figure
%         subplot(1,2,1)
        imshow(mat2gray(5*frame/max(frame(:))+squeeze(gradient(sum(newmasks,1)))));
        saveas(gcf,[pathplanes{(p-1)*num_planes+d},fileplanes{(p-1)*num_planes+d}(1:end-4),'frame.jpg'])
        figure
%         subplot(1,2,2)
        indmasks=bsxfun(@times,newmasks,(1:size(newmasks,1))');
        imagesc(squeeze(sum(indmasks,1)));
        saveas(gcf,[pathplanes{(p-1)*num_planes+d},fileplanes{(p-1)*num_planes+d}(1:end-4),'masks.jpg'])
        
        F=extract_F_poly_from_auto_ROI(newmasks,video);
        Fc=zeros(size(F));
        Fs=zeros(size(F));
        for j=1:size(F,2)
            
            junk=F(:,j);
            %     junk=junk;
            
            window=round(numframes/numwindow);
            junk2=zeros(size(junk));
            for k=1:length(junk)
                cut=junk(max(1,k-window):min(numframes,k+window));
                cutsort=sort(cut);
                a=round(length(cut)*.08);
                junk2(k)=cutsort(a);
            end
            Fc(:,j)=(junk./junk2);
            maxval=max(Fc(:,j));
            Fc(:,j)=(Fc(:,j)-median(Fc(:,j)))/max((Fc(:,j)-median(Fc(:,j))));
            Fc(:,j)=maxval*Fc(:,j);
            
            %Fc(:,i)=(junk-junk2)-mean((junk-junk2));
        end
        
        masks=newmasks;
        fullFname=[pathplanes{(p-1)*num_planes+d} [fileplanes{(p-1)*num_planes+d}(1:end-4) '_roibyauto_F.mat']];
        save(fullFname,'F','Fc','masks','numframes','N','M','numwindow');
        
    end
end








