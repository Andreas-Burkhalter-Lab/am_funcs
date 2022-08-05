function stims=parse_stim_seq(h)
% convert stimulus sequence in header to nx2 matrix
% @syntax:
%    stims=parse_stim_seq(h)
% @pre:
%    h: header return by readheader
% @post:
%    stims: nx2 matrix (1st column are valves; 2nd are time)
%           note that if file is shorter than the sequence, the 
%           sequence will be truncated.

   if(is_robot_stim(h))
      error(['parse_stim_seq() is meaningless for robot deliveried stimuli.']);
   end
   
   % get the whole stimulus sequence
   stims=eval(['[' key2value(h.wholeheader, 'stimulus sequence') ']']);
   
   if(isempty(stims)) 
      return;
   end

   % correct stims according to field nscans
   ttWholeTime=h.nscans/h.scanrate;
   if(ttWholeTime<stims(end, 2))
      % if file size is shorter than stimulus sequence
      ttIdx=find(stims(:,2)<=ttWholeTime);
      ttIdx=ttIdx(end);
      if(stims(ttIdx, 2)==ttWholeTime)
         stimsAfterTrunc=stims(1:ttIdx, :);
      else
         stimsAfterTrunc=[stims(1:ttIdx, :); stims(ttIdx, 1) ttWholeTime];
      end
      stims=stimsAfterTrunc;
   else
      if(ttWholeTime>stims(end, 2))
         warning('file size is longer than stimulus time');
      end
   end
   
   
