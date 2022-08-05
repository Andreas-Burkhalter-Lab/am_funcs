function wp = whisparms(t,sng,panal)
     % WHISPARMS: measure parameters of individual whistles
% wp = whisparms(t,sng,panal)
% where
%   t is a vector of whistle start times
%   sng is the sparse sonogram
%   panal is a structure of analysis parameters with fields
%     filts: a tensor of filters, of size filts(nfreq,ntimes,nfilts)
% and
%   wp is an output structure of measured parameters
