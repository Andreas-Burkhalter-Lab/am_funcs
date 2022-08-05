function success=reorder_merec(ifilename, ofilename, newLabelList)
% PRE:
%    newLabelList: a cell array of strings that is new label list
   
   [iheader, ifid]=readheader(ifilename);
   ofid=fopen(ofilename, 'w');
   
   oldChannelList=iheader.channels;
   oldLabelList=cellstr(split_label(key2value(iheader.wholeheader, 'label list')));
   
   tt=unique(newLabelList);
   if(length(tt)~=length(newLabelList))
      error(['reorder_merec: new label list is not unique']);
   end
   
   indices=zeros(1, length(newLabelList))-1;
   for idx=1:length(newLabelList)
      tt=strmatch(newLabelList{idx}, oldLabelList, 'exact');
      if(length(tt)~=1)
	 if(length(tt)==0)
	    error(['new label list has label ' newLabelList{idx} ' that is not in old label list']);
	 else
	    error(['label ' newLabelList{idx} ' appears multiple times in old label list']);
	 end
      end
      indices(idx)=tt;
   end
   
   % now oldLabelList(indices) is same as newLabelList
   
   strHeader=iheader.wholeheader;
   
   % TODO: update nscans. It is fine now w/o updating it b/c readheader() takes care of this.
   
   % update channel list
   newChannelList=oldChannelList(indices);
   newChannelListStr=mat2str(newChannelList);
   newChannelListStr=newChannelListStr(2:end-1); % remove [ and ]
   [strHeader, tt]=update_value(strHeader, 'channel list', newChannelListStr);
   
   % update label list
   newLabelListStr=[];
   for idx=1:length(newLabelList)
      newLabelListStr=[newLabelListStr newLabelList{idx} '$'];
   end
   [strHeader, tt]=update_value(strHeader, 'label list', newLabelListStr);
   
   % TODO: update comment. So that later we know the data is not original raw data.
   
   % write the new header
   update_header(ofid, strHeader);
   
   % a by-product of update_header() is that the file position is moved to where data should be.

   success=1;
   
   % write re-ordered data
   nScansPerBlock=1000;
   dataFormat='*int16';
   while(1)
      block=fread(ifid, [iheader.numch nScansPerBlock], dataFormat);
      if(isempty(block))
	 break;
      end
      blockNew=block(indices, :);
      count=fwrite(ofid, blockNew, dataFormat(2:end)); % 2:end: skip *
      if(count~=numel(blockNew))
	 % error('reorder_merec: write file error');
	 success=0;
	 break;
      end
   end
   
   fclose(ifid);
   fclose(ofid);
   
