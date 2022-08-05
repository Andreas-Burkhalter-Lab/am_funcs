function ranges=segment(seq, delimiters)
% segment: similar to split_str(), this func returns index ranges 
%          who define the segments separated by delimiters.
% USAGE:   
%    ranges=segment(seq, delimiters)
% PRE: 
%    seq: sequence, a vector of numbers
%    delimiters: a vector of numbers too
% POST:
%    ranges: a Rx2 matrix
%   
   ranges=[]; found=0;
   for idx=1:length(seq)
      if(any(delimiters==seq(idx)))
	 found=0;
	 continue;
      else
	 if(found)
	    ranges(end,2)=idx;
	 else
	    ranges(end+1,:)=[idx idx];
	    found=1;
	 end
      end
   end
