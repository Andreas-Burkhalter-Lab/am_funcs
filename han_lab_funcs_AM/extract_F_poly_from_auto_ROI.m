function [Fn,nF]=extract_F_poly_from_auto_ROI(masks,chonedata,varargin)
numframes=size(chonedata,3);
refimage=sum(chonedata,3);

%extract <F> in ROIs
numcolors=size(masks,1);
F=zeros(numframes,numcolors);
nF=F;
%  Fr=zeros(numframes,numcolors);

if nargin==2
    for j=1:numcolors
        regionmask=squeeze(masks(j,:,:));
        parfor i=1:numframes
            currframe=squeeze(chonedata(:,:,i));
            %         regionpixels=find(regionmask==1);
            %         cutref=refimage(regionpixels);
            cut=currframe(regionmask==1);
            F(i,j)=mean(cut);
        end
    end
    Fn=F;
else
    for j=1:numcolors
        regionmask=squeeze(masks(j,:,:));
        nmask=squeeze(varargin{1}(j,:,:));
        parfor i=1:numframes
            currframe=squeeze(chonedata(:,:,i));
            %         regionpixels=find(regionmask==1);
            %         cutref=refimage(regionpixels);
            cut=currframe(regionmask==1);
            ncut=currframe(nmask==1);            
            F(i,j)=mean(cut);
            nF(i,j)=mean(ncut);
        end
        
    end
           %Remove first PC of noise:
%             vec=pca(nF);
%             mNoise=vec(:,1)'*bsxfun(@minus,nF,mean(nF))';
%             scale=mean(nF)/max(nF(:));
%             Fn=F-(scale'*mNoise)';
            
            Fn=bsxfun(@minus,F,bsxfun(@minus,nF,mean(nF)));
end