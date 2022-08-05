function [lookup, labels]=imagine_stim_lookup(header)
% [lookup, labels]=imagine_stim_lookup(header)
% pre:
%    header: the content of .imagine file
% post:
%    lookup: a vector of valve numbers.
%    labels: a cell vector s.t. labels{v} is the descriptive name for valve v.
% 
   strStimContent=key2value(header, 'stimulus file content');
   strStimContent=strrep(strStimContent, '\n', char(10)); % decode \n to LF
   lines=cellstr(split_str(strStimContent, char(10)));    % each line goes its own cell
   
   headerSize=key2value(strStimContent, 'header size');
   headerSize=sscanf(headerSize, '%d', 1);
   
   stim=[];
   for idx=headerSize+1:length(lines)
      line=strtrim(lines{idx});
      if(isempty(line))
	 continue;
      end
      
      [numbers, count, err, nextidx]=sscanf(line, '%d', 2);
      stim=[stim; numbers'];
      label=strtrim(line(nextidx:end));
      if(numbers(1)>0)
         labels{numbers(1)}=label;
      end
   end
   
   if(isempty(stim))
      lookup=[]; labels={};
      return;
   end

   if(stim(1,2)~=0)
      error('wrong stimulus format');
   end
   
   for idx = 1:length(labels)
       if isempty(labels{idx})
           labels{idx} = ['Vlv' num2str(idx)];
       end
   end
   
   duration=diff(stim(:,2));
   
   lookup=[];
   for idx=1:length(duration)
      lookup=[lookup; ones(duration(idx), 1)*stim(idx,1)];
   end
   
   lookup(end+1)=stim(end,1);
   