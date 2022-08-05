function result=break_func_handles(list)
% break function handle/arg list into cell array of handle-arg.
% SYNTAX:
%    result=break_func_handles(list)
% EG:
%    t={@sin, @gt, 12, 34, @cos}
%    result=break_func_handles(t)
%    result{:}

   n=length(list);
   isFuncH=ones(1, n);
   for idx=1:n
      isFuncH(idx)=isa(list{idx}, 'function_handle');
   end
   
   argIndices=segment(isFuncH, 1);
   handleIndices=find(isFuncH);

   result=cell(1, length(handleIndices));
   
   for resultIdx=1:length(result)
      handleIdx=handleIndices(resultIdx);
      row=find(argIndices(:,1)==handleIdx+1);
      if(isempty(row))
         result{resultIdx}={list{handleIdx}}; % w/o any arg
      else
         result{resultIdx}=list([handleIdx argIndices(row,1):argIndices(row,2)]);
      end
   end
   
