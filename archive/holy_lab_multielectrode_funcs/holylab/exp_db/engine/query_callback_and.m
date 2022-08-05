function wanted=query_callback_and(entry, conditions)
% this is a general query callback which implements AND of subquery on
% entry's field values.
% SYNTAX:
%    wanted=query_callback_and(entry, conditions)
% PRE:
%    entry: an xdb entry
%    conditions: a struct
%       its field names are subset of entry's field names;
%       each field value is either a function handle, or
%          a cell with value {funcHandle1, extraArg11, extraArg12, ... funcHandle2, ...}
% EG:
%  conditions=struct(...
%    'investigator', {{@include, 'jason'}}, ...
%    'start_date', {{@after, datenum('20040501', 'yyyymmdd'), @before, datenum('20040601', 'yyyymmdd')}} ... 
%    );
% 
%  result=search_xdb({@query_callback_and, conditions})
%
% See also: QUERY_CALLBACK_BEHAVIOR

   wanted=0;
   if(~isstruct(conditions))
      return;
   end

   fields=fieldnames(conditions);
   % entryFields=fieldnames(entry);
   % if(~isempty(setdiff(fields, entryFields)))
   %    error('there''re unknown fields in the conditions');
   % end
   
   for idxField=1:length(fields)
      f=fields{idxField};
      
      if(~isfield(entry, f))
         return; % treat missing-fields as mismatch
      end
      
      operator=conditions.(f);
      if(iscell(operator))
         % a func handle/arg list
         subOps=break_func_handles(operator);
         for subOpIdx=1:length(subOps)
            subOp=subOps{subOpIdx};
            if(~subOp{1}(entry.(f), subOp{2:end}))
               return; % b/c not-match
            end
         end
      else
         % just a func handle
         if(~operator(entry.(f)))
            return;
         end
      end
   end % for, each field
   
   wanted=1;
   
   