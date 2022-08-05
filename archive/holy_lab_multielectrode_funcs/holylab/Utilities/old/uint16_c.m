function u=uint16_c(s)
% int16_c: cast s in C's way
% usage:
%    u=uint16_c(s)
% pre:
%    s: an int16 array (s: signed)
% post:
%    u: an uint16 array
% note:
%    1. c's way: re-interpret the 16 bits as unsigned
%    2. matlab's way: map out of range values to either 0 or 2^16-1, 
%                     depending on who is closer.
% NOTE: Use typecast instead
   
   if(~isa(s,'int16')) 
      error('int16_to_uint16: input arg is not int16'); 
   end % if, the input is not 16 bit
   
   s=double(s); % todo: no need for matlab 7+
   
   indices=find(s<0);
   s(indices)=bitcmp(-s(indices),16)+1;
   u=uint16(s);
   
