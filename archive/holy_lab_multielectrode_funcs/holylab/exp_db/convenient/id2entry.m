function result=id2entry(id)
% get an entry from its id
% SYNTAX:
%   result=id2entry(id)
% PRE:
%   id: a string, the id of the entry
% POST:
%   result: a struct, the entry

   r=search_xdb(@(entry) isequal(entry.id, id));
   if(isempty(r))
      result=[];
   else
      if(length(r)>1)
         error('database error: more than one entry have same id');
      else
         result=r{1};
      end
   end
