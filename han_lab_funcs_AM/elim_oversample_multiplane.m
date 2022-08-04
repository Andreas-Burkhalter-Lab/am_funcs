% Code to eliminate rois that were sampled in multiple planes
% find centroids <100pixels away on diff planes
clear all
planes=input('number of plane: ');
paths{planes}=0;
names{planes}=0;
Fs=4;
for p=1:planes
    disp(['Plane ', num2str(p)]);
    [names{p},paths{p}]=uigetfile('*.mat','pick your files');
end
constrain = @(sig) (sig-min(sig(:)))/max(sig-min(sig(:)));

mask{planes}=0;
dFFp{planes}=0;
p_centroids{planes}=0;
col_centroids=[];
for p=1:planes
    load([paths{p},names{p}],'dFF','F','dFF','nF','masks','M','N','Fc','frame');
    figure('units','normalized', 'Position', [.01 .05 .8 .8]);
lose=[];
for i=1:size(dFF,2)
    %Show current cell in red, 9 others in black
    subplot(1,3,1)
    plot(bsxfun(@plus,dFF,1:size(dFF,2)),'color','k')
    hold on
    plot(dFF(:,i)+i,'color','r')
    ylim([(max(1,i-10)-.5) (min(size(dFF,2),i+10)+.5)])
    %Show current cell and display skewness, peak to peak mag, and peak to
    %rms value
    subplot(1,3,2)
    plot(dFF(:,i))
    hold on
    title(['Cell ',num2str(i), ' Sk=',num2str(skewness(dFF(:,i))),...
        ' P2P=',num2str(peak2peak(dFF(:,i))),' P2R=',num2str(peak2rms(dFF(:,i)))])
    %Show 60 seconds of cell activity
    subplot(1,3,3)
        hold on
    plot(dFF(60:120*Fs,i));
%     title(['Raw F: ', num2str(mean(rawF(:,i))/mean(mean(nF)))])
%     plot(lped(60:120*Fs)+1)
    xlim([0 60*Fs])
%     [c1,d1]=ginput(1);
    keep=input('Keep? ');
    clf
    if ~keep
        lose(end+1)=i;
    end
end 
    F(:,lose)=[];
    dFF(:,lose)=[];
    nF(:,lose)=[];
    masks(lose,:,:)=[];
    Fc(:,lose)=[];
    save([paths{p},names{p}],'F','dFF','nF','masks','Fc','-append');
    mask{p}=masks;
    dFFp{p}=dFF;
    Fraw{p}=F;
    frames{p}=frame;
%     Fc{p}=Fc;
    p_centroids{p}=find_centroids(mask{p});
    col_centroids=[col_centroids; find_centroids(mask{p})];
end

%%
elim=cell(planes,1);
for p=1:planes
    for p2=p:planes
        remove=[];
        likely_overlap=[];
        cor_overlap=[];
        %Only check if both planes have cells
        if ~isempty(p_centroids{p}) && ~isempty(p_centroids{p2})
        distance=pdist2(p_centroids{p},p_centroids{p2});
        %When different planes find any centroids with 50 pixels
        if p2~=p
            possible_overlaps=find(distance<50);
        %On same plane find centroids >1 and <80 pixels away
        else
            possible_overlaps=intersect(find(distance<80),find(distance>1));
        end
        
        [cell_plane1,cell_plane2]=ind2sub(size(distance),possible_overlaps);
        corrs=zeros(size(cell_plane1,1));
        if ~isempty(cell_plane1)
        for i=1:size(cell_plane1,1)
            if p2~=p
%                 for j=1:size(cell_plane2,1)
                    % Trim traces
                    if length(dFFp{p})>length(dFFp{p2})
                       dFFp{p}=dFFp{p}(1:length(dFFp{p2}),:); 
                    else
                        dFFp{p2}=dFFp{p2}(1:length(dFFp{p}),:); 
                    end
                    rval=corrcoef(dFFp{p}(:,cell_plane1(i)),dFFp{p2}(:,cell_plane2(i)));
%                     if rval(2,1)>.9 && p==p2
%                         likely_overlap=[likely_overlap;cell_plane1(i),cell_plane2(j)];
                    if rval(2,1)>.2
                        disp(['Possibly : Plane ',num2str(p),' cell ', num2str(cell_plane1(i)),...
                            ' and plane' ,num2str(p2),' cell ', num2str(cell_plane2(i)), ' Loop Loc ', num2str(i)])
                        likely_overlap=[likely_overlap;cell_plane1(i),cell_plane2(i)];
                        cor_overlap=[cor_overlap;rval(2,1)];
                    end
%                 end
            else  
                    rval=corrcoef(dFFp{p}(:,cell_plane1(i)),dFFp{p2}(:,cell_plane2(i)));
                    if rval(2,1)>.6
                        likely_overlap=[likely_overlap;cell_plane1(i),cell_plane2(i)];
                        cor_overlap=[cor_overlap;rval(2,1)];
%                     elseif rval(2,1)>.8
%                         likely_overlap=[likely_overlap;cell_plane1(i),cell_plane2(j)];
                    end
            end
        end
        end
        if ~isempty(likely_overlap)
            figure('units','normalized', 'Position', [.01 .05 .98 .87]);
            subplot(2,2,1)
            imagesc(squeeze(sum(mask{p}(unique(likely_overlap(:,1)),:,:))))
            title(['Plane ', num2str(p)])
            subplot(2,2,2)
            imagesc(squeeze(sum(mask{p2}(unique(likely_overlap(:,2)),:,:))))
            title(['Plane ', num2str(p2)])
            subplot(2,2,[3,4])
            hold on
            plot(bsxfun(@plus,dFFp{p}(:,likely_overlap(:,1)),1:size(likely_overlap,1)))
            plot(bsxfun(@plus,dFFp{p2}(:,likely_overlap(:,2)),1:size(likely_overlap,1)))
            for ii=1:size(likely_overlap,1)
                figure('units','normalized', 'Position', [.01 .05 .98 .87]);
                subplot(2,3,1)
                imagesc(squeeze((mask{p}(likely_overlap(ii,1),:,:)))+squeeze((mask{p2}(likely_overlap(ii,2),:,:))))
                title(['Plane ', num2str(p), ' Cell ', num2str(likely_overlap(ii,1))...
                    ' ','Plane ', num2str(p2), ' Cell ', num2str(likely_overlap(ii,2)) ])
                
                m1=squeeze(mask{p}(likely_overlap(ii,1),:,:));
                m2=squeeze(mask{p2}(likely_overlap(ii,2),:,:));
                f1=frames{p};
                f2=frames{p2};
                
                [gx1,gy1]=gradient(double(m1));
                [gx2,gy2]=gradient(double(m2));
                
                pic1=zeros(size(m1,1),size(m1,2));
                pic2=pic1;
                pic1(:,:,1)=mat2gray(f1);
                pic1(:,:,2)=mat2gray(f1);
                pic1(:,:,3)=mat2gray(f1)+mat2gray(squeeze(abs(gx1)+abs(gy1)));
                
                pic2(:,:,1)=mat2gray(f2);
                pic2(:,:,2)=mat2gray(f2);
                pic2(:,:,3)=mat2gray(f2)+mat2gray(squeeze(abs(gx2)+abs(gy2)));
                
                subplot(2,3,2)
                imagesc(pic1); colormap(gray)
                title(['Plane ', num2str(p), ' Cell ', num2str(likely_overlap(ii,1)) ])
                subplot(2,3,3)
                imagesc(pic2); colormap(gray)
                title(['Plane ', num2str(p2), ' Cell ', num2str(likely_overlap(ii,2)) ])
                subplot(2,3,[4:6])
                hold on
                plot(dFFp{p}(:,likely_overlap(ii,1)))
                plot(dFFp{p2}(:,likely_overlap(ii,2)))
                title(['Correlation = ',num2str(cor_overlap(ii))]);
                legend({'Left','Right'})
                brt1=mean(Fraw{p}(:,likely_overlap(ii,1)));
                brt2=mean(Fraw{p2}(:,likely_overlap(ii,2)));
                if brt2>brt1
                    remove(ii)=ii*input(['0-Keep 1-Remove: ','(Would keep Right)']);
                    if remove(ii)
                       disp(['Removing Left Cell: Plane ',num2str(p),' Cell ',num2str(likely_overlap(ii,1))]);
                    end
                else
                    remove(ii)=ii*input(['0-Keep 1-Remove: ','(Would keep Left)']);
                    if remove(ii)
                       disp(['Removing Right Cell: Plane ',num2str(p2),' Cell ',num2str(likely_overlap(ii,2))]);
                    end
                end
                
%                 remove(ii)=ii;
                close(gcf);
            end
        end
        find(remove>0)
        if ~isempty(remove)
            for r=find(remove>0)
                if mean(Fraw{p}(:,likely_overlap(r,1)))>mean(Fraw{p2}(:,likely_overlap(r,2)))
                    elim{p2}=[elim{p2}, likely_overlap(r,2)];
                else
                    elim{p}=[elim{p}, likely_overlap(r,1)];
                end
%                 if remove(r)==2
%                     elim{p2}=[elim{p2}, likely_overlap(r,2)];
%                 else
%                     elim{p}=[elim{p}, likely_overlap(r,1)];
%                 end
            end
        end
        end
    end
end
%%
for p=1:planes
    load([paths{p},names{p}],'F','dFF','nF','masks','Fc');
    elim{p}=unique(elim{p});
    disp(['Plane ',num2str(p),' Eliminating Cells ', num2str(elim{p})]);
    F(:,elim{p})=[];
    Fc(:,elim{p})=[];
    dFF(:,elim{p})=[];
    masks(elim{p},:,:)=[];
    nF(:,elim{p})=[];
    notUsed=elim{p};
    save([paths{p},names{p}],'F','dFF','nF','masks','notUsed','Fc','-append');
end
%%
maskc=zeros(size(masks,2),size(masks,3),planes);
for p=1:planes
    load([paths{p},names{p}],'F','dFF','masks');
    if ~isempty(masks)
    if size(masks,1)>1
    maskc(:,:,p)=squeeze(sum(masks,1));
    else maskc(:,:,p)=squeeze(masks);
    end
    end
end
save([paths{p},names{p}],'maskc','-append');
