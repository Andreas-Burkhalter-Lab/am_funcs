function result=detect(video, roifile, nMice, isShowResult, timeRange, weightedMean, textureFilter)

rois=parse_rois(roifile);

result=cell(1, length(rois));

id=ffmpeg('open', video);
ffmpeg('decode', id);

duration=ffmpeg('getDuration', id);
fps=90000/3003;
nFrames=(duration-1)*fps; % ignore the data from the last sec
if(isempty(timeRange))
   frameIndexFrom=1;
   frameIndexTo=nFrames;
else
   frameIndexFrom=round(max(1, timeRange(1)*fps));
   if(timeRange(2)<=0)
      timeRange(2)=duration;
   end
   frameIndexTo=round(min(nFrames, timeRange(2)*fps));
   if(timeRange(1)>0)
      status=ffmpeg('seekPerfectFrame', id, timeRange(1));
      if(status==0)
         error('error to seek video to specified time range');
      end
   end
end

isShowInterResult=0; % for debug only
previousMasks=[];

for idx=frameIndexFrom:frameIndexTo
   frame=ffmpeg('getFrame', id); 
   frameTime=ffmpeg('getCurTime', id);
   if isShowInterResult
      figure; image(frame); title('full frame')
   end

   for roiIdx=1:length(rois)
      roi=rois(roiIdx);
      if(frameTime<roi.from || frameTime>roi.to || roi.active==0)
         continue;
      end
      % range are get from 1st and 3rd pts
      rowRange=sort([roi.geo(2) roi.geo(6)]);
      colRange=sort([roi.geo(1) roi.geo(5)]);
      % only process the roi
      subframe = frame(rowRange(1):rowRange(2),colRange(1):colRange(2),:);
      if isShowInterResult
         figure; image(subframe); axis image; title('subframe')
      end

      masks=mask_objs(subframe, nMice, isShowInterResult, textureFilter, 'g', 0.1);
      if(length(masks)<nMice)
         if(isempty(previousMasks))
            error('don''t know what to do');
         else
            masks=previousMasks;
         end
      end
      previousMasks=masks;
      
      % now masks is what we want, and we can mask out the mouse images
      % use the grayscale image so we can correlate objs across frames by sorting intensity
      % the idea is: even the abs. intensities change, but the intensity order
      %              doesn't change
      image_gray=color2gray(subframe);

      intensity=NaN(1,nMice); % avg intensity
      centers=NaN(nMice, 2);
      mm=NaN(nMice, 3);
      for mouseIndex=1:nMice
         mask=masks{mouseIndex};
         nPxl=sum(mask(:));
         intensity(mouseIndex)=sum(image_gray(mask))/nPxl; % or: sum(sum(image_gray.*mask))/nPxl;
         [rows, cols]=find(mask);
         if(weightedMean)
            centers(mouseIndex,:)=[sum(rows.*image_gray(mask))/(intensity(mouseIndex)*nPxl) ...
                                   sum(cols.*image_gray(mask))/(intensity(mouseIndex)*nPxl)];
         else
            centers(mouseIndex,:)=[mean(rows) mean(cols)];
         end
         % todo: maybe we should use each pixel's intensity as weight?
         mm(mouseIndex, 1)=sum((cols-centers(mouseIndex,2)).^2)/nPxl; % Mxx
         mm(mouseIndex, 2)=sum((rows-centers(mouseIndex,1)).^2)/nPxl; % Myy
         mm(mouseIndex, 3)=sum((cols-centers(mouseIndex,2)).*(rows-centers(mouseIndex,1)))/nPxl; % Mxy
      end
      [tt, indices]=sort(intensity, 'descend'); % masks{indices(1)} is the most intensive one
      intensity=intensity(indices);
      centers=centers(indices, :);
      mm=mm(indices, :);

      % now show the image
      if(isShowResult)
         figure; image(subframe); axis image; title(['subframe w/ center-of-mass @'  num2str(frameTime)])
         hold on
         markers={'*', 's', '^'};
         for i=1:nMice
            % index=indices(i);
            plot(centers(i,2), centers(i, 1), 'marker', markers{i});
         end
      end
      
      % now save the result in the format (assume nMice is 3):
      %  time I1 I2 I3 x1 x2 x3 y1 y2 y3 mxx1 mxx2 mxx3 myy1 myy2 myy3 mxy1 mxy2 mxy3
      result{roiIdx}(end+1,:)=[frameTime intensity centers(:)' mm(:)']; 
   end % for, each roi
   
   % for i=1:30, ffmpeg('decode', id); end
   ffmpeg('decode', id);
end % for, each frame

ffmpeg('close', id);
