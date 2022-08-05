function stringarray=split_label(astring)
% SPLIT_STR: split a string (a label list) into string list, like what split() does in perl
%                   
% Syntax:
%    stringarray=split_label(astring)
%    
% pre:   
%    astring: the string to be splited, in which labels are separated by
%             '$' or ',' or ' '.  If separated by ' ', we need concatenate
%             comment in a pair of () to previous label. Delimiter '$' is
%             considered first, then ',', then ' '.
%   
%  
% post:
%    stringarray: a matlab char array whose rows are splitted labels;
%                 return empty string if astring is empty
%    
% Notes:
%    Continuous delimiters are treated as one delimeter.

   pos_dollar=strfind(astring, '$');
   pos_comma =strfind(astring, ',');
   if(~isempty(pos_dollar)) 
      stringarray=split_str(astring, '$');
   elseif(~isempty(pos_comma))
      stringarray=split_str(astring, ',');
   else
      tStringList=split_str(astring, ' ');
      tLabels={};
      j=0;
      for i=1:size(tStringList, 1)
         if(tStringList(i,1)=='(') 
            % tLabels{j}=[tLabels{j} tStringList(i,:)];
            tLabels{j}=strcat(tLabels{j}, tStringList(i,:)); 
                  % use strcat instead of [ ] to get rid of trailling zeros
         else
            j=j+1;
            tLabels{j}= tStringList(i,:);
         end
      end
     
      stringarray=char(tLabels); % convert cell array to char array
      %    stringarray=[];
      %    for i=1:j
      %       stringarray=strvcat(stringarray, tLabels{i});
      %    end
   end
