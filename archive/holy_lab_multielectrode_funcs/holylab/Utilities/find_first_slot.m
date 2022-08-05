function slot=find_first_slot(seq)
% find the first empty slot in seq
% USAGE:
%    slot=find_first_slot(seq)
% PRE:
%    seq: a vector of natural intergers
% POST:
%    slot: smallest natural integer not in seq
   t=setdiff(1:length(seq), seq);
   if(isempty(t))
      slot=length(seq)+1;
   else
      slot=t(1);
   end
