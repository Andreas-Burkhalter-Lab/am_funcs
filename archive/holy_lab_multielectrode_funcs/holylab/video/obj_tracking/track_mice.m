function track_mice
% this is a user-friendly wrapper around detect()   
    
   videos = UIGetFiles('*.mod', 'please select video files');
   if(isempty(videos))
      return;
   end
   
   % roifile will be calc automatically
   
   prompt={'Enter #mice:'};
   name='#mice?';
   numlines=1;
   defaultanswer={'3'};
 
   answer=inputdlg(prompt,name,numlines,defaultanswer);
   nMice=str2num(answer{1});
   
   usrResponse=questdlg('Do you want to see the result for each frame?', ...
      'track mice', ...
      'yes', 'no', 'yes');
   isShowResult=isequal(usrResponse, 'yes');
   
   prompt={'Enter time range (use 0 0 for the file length):'};
   name='time range?';
   numlines=1;
   defaultanswer={'0 0'};
 
   answer=inputdlg(prompt,name,numlines,defaultanswer);
   timeRange=str2num(answer{1});
   
   [filename, pathname] = uiputfile('*.track_mice', 'Please pick a file to save the tracking result');
   if isequal(filename,0) || isequal(pathname,0)
       disp('User pressed cancel')
   end
   resultfile= fullfile(pathname, filename);
    
   for fileIdx=1:length(videos)
      trackResult(fileIdx).video=videos{fileIdx};
      trackResult(fileIdx).roifile=[trackResult(fileIdx).video '.rois'];
      trackResult(fileIdx).result=detect(trackResult(fileIdx).video, trackResult(fileIdx).roifile, ...
         nMice, isShowResult, timeRange, 1, @stdfilt);
   end
   
   save(resultfile, 'trackResult');
   