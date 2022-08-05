function matched=query_callback_behavior(entry)
% QUERY_CALLBACK_BEHAVIOR: a sample function for searching the expt. database
%
% Query callbacks have the following syntax:
%   matched = query_callback_behavior(entry)
% where
%   entry is an xdb entry
% and
%   matched is boolean (true if the entry matches, false otherwise).
%
% For this sample query callback function, it returns true if the
% experiment used male urine as one of stimuli, unless there were any
% stimulus mixtures used in the experiment.
%
% See also: QUERY_CALLBACK_AND.

   matched=false;
   
   stimuli=entry.stimulus;
   
   if(isempty(stimuli)) return; end
   
   for idxStim=1:length(stimuli)
      stim=stimuli{idxStim};
      if(length(stim)>1)
        % length > 1 indicates a mixture
         return; 
      end % mixture
      
      if(stim.identity=='M') matched=true; end
   end
