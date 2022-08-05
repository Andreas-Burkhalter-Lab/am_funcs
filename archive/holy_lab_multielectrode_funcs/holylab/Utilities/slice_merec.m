function success=slice_merec(ofilename, inputInfo)
% copy parts of some .merec files and form a new .merec file
% SYNTAX:
%    success=slice_merec(ofilename, inputInfo)
% PRE:
%    ofilename: output file name
%    inputInfo: a struct array which has two fields:
%       .filename: the filename to get data from
%       .range: a 1x2 matrix, which specifies time ranges (in sec)
% NOTE: 
%    It is assumed that all input files:
%       have exactly same input channel list (include channel order).
   
   [iheader, ifid]=readheader(inputInfo(1).filename);
   nOutputScans=0;
   for idx=1:length(inputInfo)
      nScansCur=diff(time2Scans(inputInfo(idx)))+1;
      nOutputScans=nOutputScans+nScansCur;
   end
   
   ofid=fopen(ofilename, 'w');
      
   strHeader=iheader.wholeheader;
   
   % update nscans
   [strHeader, tt]=update_value(strHeader, 'nscans', num2str(nOutputScans));
   
   % update comment, So that later we know the data is not original raw data.
   newComment='';
   for idx=1:length(inputInfo)
      newComment=sprintf('%s%s:[%s]; ', newComment, inputInfo(idx).filename, ...
         num2str(inputInfo(idx).range));
   end
   [strHeader, tt]=update_value(strHeader, 'comment', newComment);
   
   % write the new header
   update_header(ofid, strHeader);
   
   % a by-product of update_header() is that the file position is moved to where data should be.

   fclose(ifid);
   
   success=1;
   
   % write data
   blockSize=1000;
   dataFormat='*int16';
   for idxInput=1:length(inputInfo)
      memm = merecmm(inputInfo(idxInput).filename,'tovolts',false,'contiguous',true);
      iheader=memm.header;
      channels=iheader.channels;
      scanRange=time2Scans(inputInfo(idxInput));
      nBlocks=ceil((diff(scanRange)+1)/blockSize);
      for idxBlock=1:nBlocks
         scanNumFrom=(idxBlock-1)*blockSize+scanRange(1);
         scanNumTo=min(scanNumFrom+blockSize-1, scanRange(2));
         wave= memm(channels, [scanNumFrom scanNumTo]); % each col is a scan
         
         count=fwrite(ofid, wave, dataFormat(2:end)); % 2:end: skip *
         if(count~=numel(wave))
            success=0;
            break;
         end
      end % while(1)
      if(~success)
         break;
      end
   end % for, each inputInfo
   
   fclose(ofid);
   
   
% convert time  range to scan range   
function scanRange=time2Scans(inputInfo)
   filename=inputInfo.filename;
   range=inputInfo.range;
   header=readheader(filename);
   scanRange=range*header.scanrate;

   scanRange(1)=scanRange(1)+1; % +1: matlab is 1-base numbering.
   if(scanRange(1)>header.nscans)
      error('tried to slice beyond file end');
   end
   if(scanRange(2)<0 || scanRange(2)>header.nscans)
      scanRange(2)=header.nscans;
   end
   
   
   