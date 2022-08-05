videofile='/usr/lab/video/MC13_base_5.12.10.mov';
id=ffmpeg('open', videofile);

status=ffmpeg('seekPerfectFrame', id, 32); % seek to 32 sec
ffmpeg('decode', id); % read packets until a frame is decodable then decode them into a frame

clear allframes
for i=1:10
   frame=ffmpeg('getFrame', id); % read the decoded frame
   frameTime=ffmpeg('getCurTime', id);
   % subframe = frame(rowRange(1):rowRange(2),colRange(1):colRange(2),:);

   allframes(:,:,:,i)=frame;
   ffmpeg('decode', id);
end
ffmpeg('close', id);

figure; image(allframes(:,:,:,1)); axis image; title('frame #1')


