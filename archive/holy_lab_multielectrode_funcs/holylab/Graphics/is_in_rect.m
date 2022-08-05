function result=is_in_rect(pt, rect, include_border)
% Syntax:
%   result = is_in_rect(pt,rect)
% or
%   result = is_in_rect(pt,rect,include_border)
%
% rect: in the form of (x1, y1, x2, y2)
% pt: in the form of (x,y)
%   pt can be an n-by-2 matrix, in which case result is a vector
   if(nargin==2) include_border=1; end
   product1=(pt(:,1)-rect(1)).*(pt(:,1)-rect(3));
   product2=(pt(:,2)-rect(2)).*(pt(:,2)-rect(4));
   if(include_border)
      result = product1<=0 & product2<=0;
   else
      result = product1<0 & product2<0;
   end
   
