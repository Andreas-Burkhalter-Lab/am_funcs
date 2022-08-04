function autoROI(pathplane,fileplane,pathref,planeref)
[filename,pathname]=uigetfile('Motion corrected plane to add rois');
load([pathname,filename]);

[filename,pathname]=uigetfile('Select Reference Frame');
maskref=load([pathname,filename],'masks');
ref_centroids=find_centroids(maskref.masks);
numref=size(maskref,1);
newmasks=zeros(size(maskref);
frame=mean(video,3);
rois=bwlabel(frame>(mean(frame(:)+2*std(frame(:)))));
rinfo=regionprop(rois);
num_rois=length(rinfo);
centroids=cat(1,rinfo.Centroid);
area=cat(1,rinfo.Area);
sizeable=find(area>150);
distmat=pdist2(ref_centroids,centroids(area>150));
topleft=find(distmat(1,:)==min(distmat(1,:)));

ny=centroids(sizeeable(topleft),1);
nx=centroids(sizeeable(topleft),2);

distx=nx-ref_centroids(1,2);
disty=ny-ref_centroids(1,1);

for i=1:numref
    tempmask=zeros(size(frame));
    m=squeeze(maskref(i,:,:));
    [r,c]=find(m>0);
    r=r+ny;
    c=c+nx;
    
    locs=sub2ind(size(frame),r,c);
    tempmask(locs)=1;
    newmasks(i,:,:)=tempmask;
end
numframes=length(chone);
N=size(chone,2);
M=size(chone,3);
F=extract_F_poly_from_auto_ROI(newmasks,video);
    Fc=zeros(size(F));
    Fs=zeros(size(F));
    for j=1:size(F,2)

        junk=F(:,j);
        junk=junk;

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
       
    
 fullFname=[path [filename(1:end-4) '_roibyhand_F.mat']];
 save(fullFname,'F','Fc','colors','masks','numframes','N','M','numwindow');










