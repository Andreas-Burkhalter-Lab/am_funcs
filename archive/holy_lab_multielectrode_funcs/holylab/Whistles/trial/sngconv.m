function convo= sngconv(sng,filt)
     % SNGCONV: convolve a sparse sonogram with a detection filter
% convo = sngconv(sng,filt)
% where
%   sng is a sparse sonogram
%   filt is the filter use
% and
%   convo is the result of the convolution
%
% See also: BUILDBAR, SPARSESNG, SPSNGPLOT.
convo = zeros(1,size(sng,2)-size(filt,2)+1);
if ~isreal(sng)
    sng = abs(sng);
end
for i = 1:length(convo)
     convo(i) = sum(sum(filt .* sng(:,i:i+size(filt,2)-1)));
end


