function [waveform] = createwaveform(nscans,tsnip,spike)
% CREATEWAVE: Creates raw waveforms
%
% This function creates raw waveforms with spikes. The spikes are
% interpolated if they occur at non-integer times.
%
% Syntax:
% createwaveform(nscans,tsnip,spike)
% where
%   nscans is a scalar specifying the number of scans in the waveform
%   tsnip is a cell array of vectors containing the times (scan numbers)
%     when the spikes occur.  All the times in a particular vector 
%     correspond to the same spike shape.
%   spike is a cell array of vectors containing points which specify
%     the shapes of the spikes.


%  nscans = 1000;
%  tsnip = {[100,200.4 300.6 400],[500 700.5 900]};
%  spike = {[500 1000 500 -200],[1000 -200]};

if (~iscell(tsnip))
    tsnip = {tsnip};
end

if (~iscell(spike))
    spike = {spike};
end

if (length(spike) ~= length(tsnip))
    errordlg('The input variables tsnip and spike must contain the same number of vectors','Input Error');
    return 
end    
    
waveform = zeros(1,nscans);

snippet = {};
newsnippet = {};
x = {};

for j = 1:length(tsnip)
    for i = 1:length(tsnip{j})
        snippet{i} = [0 spike{j}];
        snippet{i}(end + 1) = 0;
        x{i} = round(tsnip{j}(i)) - tsnip{j}(i);
        newsnippet{i} = 0.5*x{i}*(x{i}-1)*snippet{i}(1:end-2) - (x{i}-1)*(x{i}+1)*snippet{i}(2:end-1) + 0.5*x{i}*(x{i}+1)*snippet{i}(3:end);
        for k = 0:length(newsnippet{i})-1
            waveform(round(tsnip{j}(i)) + k) = waveform(round(tsnip{j}(i)) + k) +  newsnippet{i}(k+1);
        end
    end
end

plot(waveform)

