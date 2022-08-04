%% align_running  - fixes sampling differences between Fcs and abf
% input (pnum) or ({paths},{names},{pnum}
function [Fall,y2,f2,r2,cell_per_plane,mask,Fraw,rewsmall2,spksout]=...
    align_running(varargin)
if nargin==1
    pnum=varargin{1};
    for p=1:pnum
        [names{p},paths{p}]=uigetfile('*.mat');
    end
else
    paths=varargin{1};
    names=varargin{2};
    pnum=varargin{3};
end
load([paths{1},names{1}]);

if exist('Fca','var')
    Fall=Fca;
else Fall=dFF;
end
if ~exist('masks','var')
    masks=0;
end
if ~exist('F','var')
    F=Fca;
end
if length(ybinned) ~= size(Fall,1)
    rewratio=length(ybinned)/(size(Fall,1)*pnum);
    forwardvel=fix_trace(forwardvel);
    rotationvel=fix_trace(rotationvel);
    for jj=1:(length(Fall)*pnum)
        if (jj*rewratio)<length(ybinned)
            rotationvel(jj)=max(rotationvel(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
            forwardvel(jj)=max(forwardvel(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
            ybinned(jj)=max(ybinned(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
        else
            rotationvel(jj)=max(rotationvel(round(((jj-1)*rewratio+1)):end));
            forwardvel(jj)=max(forwardvel(round(((jj-1)*rewratio+1)):end));
            ybinned(jj)=max(ybinned(round(((jj-1)*rewratio+1)):end));
        end
    end
    rotationvel((jj+1):end)=[];
    forwardvel((jj+1):end)=[];
    forwardvel=smooth(forwardvel);
    rotationvel=smooth(rotationvel);
    ybinned((jj+1):end)=[];
end
Fall=[];
Fraw=[];
spksout=[];
if pnum==1
    Fall=dFF;
    Fraw=F;
    y2=repmat(ybinned,1,size(dFF,2));
    f2=repmat(forwardvel,1,size(dFF,2));
    r2=repmat(rotationvel,1,size(dFF,2));
    cell_per_plane=size(dFF,2);
    mask{1}=masks;
    if exist('spks','var')
        spksout=spks;
    end
end
rewratio=length(rewards)/length(dFF);
rewsmall2=zeros(length(dFF),1);
num_cells=size(dFF,2);
for jj=1:length(dFF)
    rewsmall2(jj)=max(rewards(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
end
ysave=ybinned;
fsave=forwardvel;
rsave=rotationvel;
if pnum>1
    cell_per_plane(4)=0;
    mask{4}=0;
    for p=1:pnum
        
        load([paths{p},names{p}]);
    if ~isempty(dFF)
        if exist('dFF','var')
            Ft=dFF;
        elseif exist('Fca','var')
            Ft=Fca;    F=Fca;
        end
        if ~exist('masks','var')
            masks=0;
        end
        
        mask{p}=masks;
        cell_per_plane(p)=size(Ft,2);
        if ~isempty(Fall)
            if length(Fall)<length(Ft)
                Fall=[Fall,Ft(1:length(Fall),:)];
                Fraw=[Fraw,F(1:length(Fall),:)];
            else
                Fall=[Fall(1:length(Ft),:), Ft];
                Fraw=[Fraw(1:length(Ft),:), F];
            end
            if exist('spks','var')
%                 spksout=[spksout, spks];
            if length(spksout)<length(spks)
                spksout=[spksout,spks(1:length(spksout),:)];
            else
                spksout=[spksout(1:length(spks),:), spks];
            end
            end
        else
            Fall=Ft;
            Fraw=F;
            
            if exist('spks','var')
                spksout=spks;
            end
        end
    else
        cell_per_plane(p)=0;
    end
    end
    
    ybinned=reshape(ysave(1:(end-mod(length(ysave),pnum)))',pnum,[])';
    forwardvel=reshape(fsave(1:end-mod(length(fsave),pnum))',pnum,[])';
    rotationvel=reshape(rsave(1:end-mod(length(rsave),pnum))',pnum,[])';
    
    y2=zeros(size(Fall));
    f2=zeros(size(Fall));
    r2=zeros(size(Fall));
    for p=1:pnum
        if p==1
            y2(:,1:cell_per_plane(1))=repmat(ybinned(1:length(Fall),1),1,cell_per_plane(1));
            r2(:,1:cell_per_plane(1))=repmat(rotationvel(1:length(Fall),1),1,cell_per_plane(1));
            f2(:,1:cell_per_plane(1))=repmat(forwardvel(1:length(Fall),1),1,cell_per_plane(1));
        else
            y2(:,(sum(cell_per_plane(1:(p-1)))+1):sum(cell_per_plane(1:(p-1)))+cell_per_plane(p))=...
                repmat(ybinned(1:length(Fall),p),1,cell_per_plane(p));
            r2(:,(sum(cell_per_plane(1:(p-1)))+1):((sum(cell_per_plane(1:(p-1))))+cell_per_plane(p)))=...
                repmat(rotationvel(1:length(Fall),p),1,cell_per_plane(p));
            f2(:,(sum(cell_per_plane(1:(p-1)))+1):((sum(cell_per_plane(1:(p-1))))+cell_per_plane(p)))=...
                repmat(forwardvel(1:length(Fall),p),1,cell_per_plane(p));
        end
        
    end
    
end
end

function for1=fix_trace(forwardvel)
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