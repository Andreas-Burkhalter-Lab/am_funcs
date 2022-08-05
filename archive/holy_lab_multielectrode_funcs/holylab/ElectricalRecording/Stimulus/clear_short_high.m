function validTrans=clear_short_high(trans, thres)
% clear_short_high: get rid of short high transitions
% pre:
%    trans: mx2 matrix. first col is zeros and ones; 2nd col is time in #scans
%    thres: the time threshold in #scans
% post:
%    validTrans: the result after removing the short highs
% 
% note: first row must begin w/ low
   if(size(trans,1)==0)
      validTrans=[];
      return ;
   end
   
   % thres=h.scanrate/1000; % 1ms
   
   nTransToCheck=floor(size(trans,1)/2)-1;
   for idxHigh=1:nTransToCheck
      if(trans(idxHigh*2+1,2)-trans(idxHigh*2,2)<thres)
         trans(idxHigh*2,1)=0; % treat it as low
      end
   end

   % remove continuous low or continuous high
   validTrans=trans(1,:);
   for idxRow=2:size(trans,1)
      if(trans(idxRow,1)~=validTrans(end,1))
         validTrans=[validTrans; trans(idxRow,:)];
      end
   end
