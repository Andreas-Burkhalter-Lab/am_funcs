function header = read_merec_header(file)
% READ_MEREC_HEADER: read Merec header into a matlab structure
%                    
% Syntax:
%    header = read_merec_header(file)
%    
% pre:   
%    file: file name
%    
% post:
%    header: a matlab structure holds all Merec header info
%  
% Deprecated/Obsolete:
%    read_merec_header_cpp(fd)                        
% 
% Notes:
%    this func is intended to replace the c++ version function  read_merec_header_cpp.mexglx
%    
  
   % assume caller is sure the file is Merec data file, otherwise here we
   % should test if the file is Merec data file by checking magic number at
   % the beginning of the file

   tstrH=read_ascii_header(file);
   header.headersize =str2num(key2value(tstrH, 'header size'));
   header.nscans     =str2num(key2value(tstrH, 'nscans'));
   tstrChList=key2value(tstrH, 'channel list');
   tChList= split_str(tstrChList,' ');
%    tChList= char(strsplit(tstrChList));       %% AM switched split_str to built in functions char(strsplit()) 1/26/15
   while all(isspace(tChList(end,:)))       %% AM added to remove trailing blank rows 1/26/15
       tChList(end,:) = [];
   end
   header.numch      =size(tChList, 1); %get row #, i.e. #strings
   
   % header.channels   =str2num(tChList); %unfortunately, this doesn't
   %                                       work very well when channels
   %                                       have diff num of digits
   for i=1:size(tChList,1) header.channels(i)=str2num(tChList(i,:)); end  %this works
   % header.channels=str2double(tChList); %no luck, still doesn't work
   
   header.scanrate   =str2num(key2value(tstrH, 'scan rate'));
   
   minSample=str2num(key2value(tstrH, 'min sample'));
   maxSample=str2num(key2value(tstrH, 'max sample'));
   minInput=str2num(key2value(tstrH, 'min input'));
   maxInput=str2num(key2value(tstrH, 'max input'));
   header.scalemult  =(maxInput-minInput)/(maxSample-minSample);
   header.scaleoff   =minInput;
   header.voltageMin =minInput;
   header.voltageMax =maxInput;
   header.minSample = minSample;
   header.maxSample = maxSample;
   
   datetime=key2value(tstrH, 'datetime');
   [s,f,t]=regexp(datetime, '([^, ]*)[, ](.*)', 'once'); % find the first ,
                                                         % or space
   if(datenum(version('-date'))>=732704)
     t={t};
   end
   header.date=datetime(t{1}(1,1):t{1}(1,2));
   header.time=datetime(t{1}(2,1):t{1}(2,2));
%  datetime=split(datetime, ', '); %this has problem when spaces appear
%                                  % several times
%  header.date=datetime(1,:);
%  header.time=datetime(2,:);
   
   
   header.usrhdr=key2value(tstrH, 'comment');
   tEndian=str2num(key2value(tstrH, 'is little endian'));
   if(tEndian>0) 
     header.endian='l';
   else 
     header.endian='b';
   end

   %  now all traditional fields are filled.
   %  now are new fields:
   header.wholeheader=tstrH;
   
   tstrTimeMarks=key2value(tstrH, 'time mark');
   if(~isempty(tstrTimeMarks))
      tt=split_str(tstrTimeMarks, ';');
      tt=cellstr(tt);
      for idxTimeMark=1:length(tt)
         tStr=tt{idxTimeMark};
         tPos=strfind(tStr, ':');
         tPos=tPos(1);
         tTimeMark.time(idxTimeMark)=str2num(tStr(1:tPos-1));
         tTimeMark.comment{idxTimeMark}=tStr(tPos+1:end); % todo: decode \n and \s
      end
      header.timemark=tTimeMark;
   else
      header.timemark=[];
   end
   
