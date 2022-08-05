function [tacout,indxout] = autocorrspikecat(t,tmax,nbins)
% autocorrspikecat: compute spike autocorrelations for vectors in a cell array
% Calling syntax is just like autocorrspike, except t may be a cell array
% of spike time vectors
binning = 0;
if (nargin == 3)
        binning = 1;
end
if (nargout > 1 & binning)
        error('Only one output when binning');
end
if (~iscell(t))                        % Allow vector inputs for generality
        if (binning)
                tacout = autocorrspike(t,tmax,nbins);
        else
                [tacout,indxout] = autocorrspike(t,tmax);
        end
        return
end
if (binning)
        tacout = zeros(1,nbins);
        for i = 1:length(t)
                tacout = tacout + autocorrspike(t{i},tmax,nbins);
        end
else
        tacout = [];
        for i = 1:length(t)
                [tactemp,indxout{i}] = autocorrspike(t{i},tmax);
                tacout(end+1:end+length(tactemp)) = tactemp;
        end
end
