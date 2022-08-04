function [Fall,y2,f2,r2,cell_per_plane,mask,Fraw,rewsmall2,nF,x2]=align_running(paths,names,pnum)
nF=[];
load([paths{1},names{1}]);

if length(ybinned) ~= size(Fc,1)
    rewratio=length(ybinned)/(size(Fc,1)*pnum);
    forwardvel=fix_trace(forwardvel);
    rotationvel=fix_trace(rotationvel);
    [forwardvel,rotationvel]=speed_transform(forwardvel,rotationvel);
    for jj=1:(length(Fc)*pnum)
        if (jj*rewratio)<length(ybinned)
            rotationvel(jj)=max(rotationvel(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
            forwardvel(jj)=max(forwardvel(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
            xbinned(jj)=max(xbinned(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
            ybinned(jj)=max(ybinned(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
            
        else
            rotationvel(jj)=max(rotationvel(round(((jj-1)*rewratio+1)):end));
            forwardvel(jj)=max(forwardvel(round(((jj-1)*rewratio+1)):end));
            xbinned(jj)=max(xbinned(round(((jj-1)*rewratio+1)):end));
            ybinned(jj)=max(ybinned(round(((jj-1)*rewratio+1)):end));
        end
    end
    rotationvel((jj+1):end)=[];
    forwardvel((jj+1):end)=[];
    forwardvel=smooth(forwardvel);
    rotationvel=smooth(rotationvel);
    xbinned((jj+1):end)=[];
    ybinned((jj+1):end)=[];
end
Fall=[];
Fraw=[];
if pnum==1
    Fall=Fc;
    Fraw=F;
    y2=repmat(ybinned,1,size(Fc,2));
    f2=repmat(forwardvel,1,size(Fc,2));
    r2=repmat(rotationvel,1,size(Fc,2));
    x2=repmat(xbinned,1,size(Fc,2));
    cell_per_plane=size(Fc,2);
    mask{1}=masks;
end
rewratio=length(rewards)/length(Fc);
rewsmall2=zeros(length(Fc),1);
num_cells=size(Fc,2);
for jj=1:length(Fc)
    rewsmall2(jj)=max(rewards(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
end

if pnum>1
    cell_per_plane(4)=0;
    mask{4}=0;
    for p=1:pnum
        load([paths{p},names{p}],'F','Fc','masks');
        mask{p}=masks;
        cell_per_plane(p)=size(Fc,2);
        if ~isempty(Fall)
            if length(Fall)<length(Fc)
                Fall=[Fall,Fc(1:length(Fall),:)];
                Fraw=[Fraw,F(1:length(Fall),:)];
            else
                Fall=[Fall(1:length(Fc),:), Fc];
                Fraw=[Fraw(1:length(Fc),:), F(1:length(Fc),:)];
            end
        else
            Fall=Fc;
            Fraw=F;
        end
    end
    
    ybinned=reshape(ybinned',pnum,[])';
    ybinned=ybinned(1:length(Fall),:);
    forwardvel=reshape(forwardvel',pnum,[])';
    forwardvel=forwardvel(1:length(Fall),:);
    rotationvel=reshape(rotationvel',pnum,[])';
    rotationvel=rotationvel(1:length(Fall),:);
    
    y2=zeros(size(Fall));
    f2=zeros(size(Fall));
    r2=zeros(size(Fall));
    x2=zeros(size(Fall));
    for p=1:pnum
        if p==1
            y2(:,1:cell_per_plane(1))=repmat(ybinned(:,1),1,cell_per_plane(1));
            r2(:,1:cell_per_plane(1))=repmat(rotationvel(:,1),1,cell_per_plane(1));
            f2(:,1:cell_per_plane(1))=repmat(forwardvel(:,1),1,cell_per_plane(1));
            x2(:,1:cell_per_plane(1))=repmat(xbinned(:,1),1,cell_per_plane(1));
            
        else
            y2(:,(sum(cell_per_plane(1:(p-1)))+1):sum(cell_per_plane(1:(p-1)))+cell_per_plane(p))=...
                repmat(ybinned(:,p),1,cell_per_plane(p));
            r2(:,(sum(cell_per_plane(1:(p-1)))+1):((sum(cell_per_plane(1:(p-1))))+cell_per_plane(p)))=...
                repmat(rotationvel(:,p),1,cell_per_plane(p));
            f2(:,(sum(cell_per_plane(1:(p-1)))+1):((sum(cell_per_plane(1:(p-1))))+cell_per_plane(p)))=...
                repmat(forwardvel(:,p),1,cell_per_plane(p));
            x2(:,(sum(cell_per_plane(1:(p-1)))+1):sum(cell_per_plane(1:(p-1)))+cell_per_plane(p))=...
                repmat(xbinned(:,p),1,cell_per_plane(p));
        end
        
    end
    
end
end

function for1=fix_trace(forwardvel)
if size(forwardvel,2)>size(forwardvel,1)
    forwardvel=forwardvel';
end
for1=forwardvel;
d=[0;diff(forwardvel)];
dd=[0;diff(d)];
for1(abs(dd)<.05)=NaN;
for2=for1;

for2(isnan(for2))=[];
for2=smooth(for2,9);
d2=[0;diff(for2)];
for2(abs(d2)>.03)=NaN;
for3=for2;

for3(isnan(for3))=[];
d3=[0;diff(for3)];
for3(abs(d3)>.03)=NaN;
for4=for3;

for4(isnan(for4))=[];
for4=smooth(for4);
d4=[0;diff(for4)];
for5=for4;
for5(abs(d4)>.03)=NaN;
tidx=1:length(for4);
for5=real(interp1(tidx(~isnan(for5)),for5(~isnan(for5)).^11,tidx,'pchip','extrap').^(1/11))';
for3(~isnan(for3))=for5;
tidx=1:length(for3);
for3=real(interp1(tidx(~isnan(for3)),for3(~isnan(for3)).^11,tidx,'pchip','extrap').^(1/11))';
for2(~isnan(for2))=for3;
tidx=1:length(for2);
for2=real(interp1(tidx(~isnan(for2)),for2(~isnan(for2)).^11,tidx,'pchip','extrap').^(1/11))';
for1(~isnan(for1))=for2;
tidx=1:length(for1);
for1=real(interp1(tidx(~isnan(for1)),for1(~isnan(for1)).^11,tidx,'pchip','extrap').^(1/11))';
end