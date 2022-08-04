 function disp_rois(varargin)
if nargin==0
nplanes=input('Number of planes ');

for p=1:nplanes
    [namevid{p},pathsvid{p}]=uigetfile('*.mat','pick your VIDEO file');
    [nameF{p},pathsF{p}]=uigetfile('*.mat','pick your ROI file');
end
else 
    nplanes=1;
    pathsvid{1}=varargin{1};
    namevid{1}=varargin{2};
    nameF{1}=varargin{3};
    pathsF{1}=varargin{1};
    
   
end
%%
prev_sum=0;
for p=1:nplanes
   load([pathsvid{p},namevid{p}])
   if exist('chone_corr','var')
       video=chone_corr;
       clear chone_corr
   end
   load([pathsF{p},nameF{p}])
    frame=std(single(video),[],3);
    
    filt=fspecial('disk',20);
    blurred=imfilter(frame,filt,'replicate');
    frame=frame./blurred;
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    if exist('A_or','var')
        
       for ii=size(A_or,2)
           
       end
    end
    masks(masks>0)=1;
    [gx1,gy1]=gradient(sum(masks(:,:,:),1));
    pic=zeros(size(masks,2),size(masks,3),3);
    pic(:,:,1)=mat2gray(frame);
    pic(:,:,2)=mat2gray(frame);
    pic(:,:,3)=mat2gray(frame)+mat2gray(squeeze(abs(gx1)+abs(gy1)));
    %         imshow(mat2gray(100*frame/max(frame(:))+squeeze(10*gradient(sum(masks(1:(end-1),:,:),1)))));
%     pic2=insertText(pic,find_centroids(masks),prev_sum+(1:(size(masks,1))),'FontSize', 8, 'BoxColor', 'White', 'BoxOpacity', 0, 'AnchorPoint', 'Center','TextColor',[.1 .8 .1]);
%     imagesc(pic)
    cents=find_centroids(masks);
%     for ii=1:size(masks,1)
%        text(cents(ii,1),cents(ii,2),num2str(prev_sum+ii),'Color','w','HorizontalAlignment','Center')
%     end
%     saveas(gcf,[pathsvid{p},namevid{p}(1:(end-4)),'AddRoisCalled.jpg'])
%     
%         figure('units','normalized', 'Position', [.01 .05 .98 .87]);
        imagesc(pic)
    for ii=1:size(masks,1)
       text(cents(ii,1),cents(ii,2),num2str(ii),'Color','w','HorizontalAlignment','Center')
    end
    saveas(gcf,[pathsvid{p},namevid{p}(1:(end-4)),'RoisCalled.jpg'])
%     prev_sum=prev_sum+size(masks,1);
    
        figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    masks(masks>0)=1;
    [gx1,gy1]=gradient(sum(masks(:,:,:),1));
    pic=zeros(size(masks,2),size(masks,3),3);
    pic(:,:,1)=mat2gray(frame);
    pic(:,:,2)=mat2gray(frame);
    pic(:,:,3)=mat2gray(frame)+mat2gray(squeeze(abs(gx1)+abs(gy1)));
            imshow(mat2gray(100*frame/max(frame(:))+squeeze(10*gradient(sum(masks(1:(end-1),:,:),1)))));
    pic2=insertText(pic,find_centroids(masks),prev_sum+(1:(size(masks,1))),'FontSize', 8, 'BoxColor', 'White', 'BoxOpacity', 0, 'AnchorPoint', 'Center','TextColor',[.1 .8 .1]);
    imagesc(pic)
    cents=find_centroids(masks);
    for ii=1:size(masks,1)
       text(cents(ii,1),cents(ii,2),num2str(prev_sum+ii),'Color','w','HorizontalAlignment','Center')
    end
    saveas(gcf,[pathsvid{p},namevid{p}(1:(end-4)),'AddRoisCalled.jpg'])
        imagesc(pic)
    for ii=1:size(masks,1)
       text(cents(ii,1),cents(ii,2),num2str(ii),'Color','w','HorizontalAlignment','Center')
    end
    saveas(gcf,[pathsvid{p},namevid{p}(1:(end-4)),'RoisCalled.jpg'])
    prev_sum=prev_sum+size(masks,1);
    
            figure('units','normalized', 'Position', [.01 .05 .98 .87]);
            subplot(1,2,1)
  imagesc(pic)
    for ii=1:size(masks,1)
       text(cents(ii,1),cents(ii,2),num2str(ii),'Color','w','HorizontalAlignment','Center')
    end
    
    subplot(1,2,2)
    plot(bsxfun(@plus,bsxfun(@rdivide,Fc,max(Fc,2)),1:size(Fc,2)))
        saveas(gcf,[pathsvid{p},namevid{p}(1:(end-4)),'RoiActivityCalled.jpg'])

end


end