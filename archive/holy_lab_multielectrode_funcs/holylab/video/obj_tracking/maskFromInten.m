function result=maskFromInten(maskFromTexture, channel, isShowInterResult)
% calc new mask based on maskFromTexture, now use intensity as the hint
   if(nargin==2)
      isShowInterResult=0;
   end
   
   mask=maskFromTexture;
   ch=double(channel);
   masked=ch(mask);
   m=mean(masked); s=std(masked);
   
   if(m>255*0.95) % NOTE: hard coded
      result=mask; % NOTE: white area. we can also use fold=1, see below
      return;
   end
   
   fold=1;
   newMask=ch>m-fold*s & ch<m+fold*s;
   newMask_closed=imclose(newMask, strel('disk',2)); % gaussian filter may also work but this is faster and return bw
   objLabel=bwlabel(newMask_closed);
   if isShowInterResult
      figure; rgb=label2rgb(objLabel); imshow(rgb); axis image; title('labeled objs from inten')
   end
   
   nObj=max(objLabel(:));
   submasks=cell(1,nObj);
   for objIndex=1:nObj
      submasks{objIndex}=objLabel==objIndex;
   end
   
   % areas=cellfun(@(x)sum(x(:)), submasks);
   for subIndex=1:length(submasks)
      submask=submasks{subIndex}; % submask is a bad name: ideally "mask" should be fully covered by a "submask"
      % change enlargeFactor based on m
      enlargeFactor=8; % NOTE: hard coded
      if(sum(submask(:))/sum(mask(:))>enlargeFactor) 
         % must be a shadow/background
         continue;
      end
      
      comm=submask & mask;
      diff=xor(comm, mask);
      relDiff=sum(diff(:))/sum(mask(:));
      if(relDiff<0.1) % NOTE: hard coded 
         result=submask;
         if isShowInterResult
            figure; rgb=label2rgb(submask); imshow(rgb); axis image; title('picked inten mask')
         end
         return;
      end
   end
   
   % don't know what to do
   result=mask;
   