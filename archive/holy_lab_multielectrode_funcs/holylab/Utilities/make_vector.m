function result=make_vector(input, type)
% result=make_vector(input, type)
% PRE:
%    input: a matrix
%    type: 'row' or 'col'
% POST:
%    result: a row or column vector
   if(strcmp(type, 'row'))
      result=reshape(input,1,numel(input));
   else
      result=reshape(input,numel(input),1);
   end
   
   
