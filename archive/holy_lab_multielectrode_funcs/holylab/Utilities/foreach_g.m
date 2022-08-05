function varargout = foreach_g(funcHandle, varargin)
%  varargout = foreach_g(funcHandle, varargin)
%  similar to deal each_return_value=funcHandle(each_input_value)
%  eg:
%    [a, b]=foreach_g(@log10, 10, 1000)
%    [a, b, c]=foreach_g(@log10, 10, 10000)
%    [a, b, c]=foreach_g(@log10, 10)
%    ca={1, 20, 300}; cb=cell(1,3);
%    [cb{:}]=foreach_g(@num2str, ca{:})
%
% Note: in many applications it's more efficient to do the following:
%     output = cellfun(funcHandle, cellarray)
%
% See also: cellfun.

   if(nargin<2)
      error('2 or more parameters are needed');
   end
   
   for idx=1:nargout
      tt=varargin{min(idx,end)}; 
      varargout{idx}=funcHandle(tt);
   end
