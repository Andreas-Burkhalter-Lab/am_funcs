function [factor,vpk] = micconvert(digital,analog,gain)

% digital: enter bit rate. If full-scale values, then enter 'fs'
% analog: 'mV','V', or 'Pa'
% gain: position of gain knob (1 - 11)
% Setting gain = 0 outputs factor of 1.0 and no conversion is performed.
% Conversion appropriate for CM16 microphone only!!!


if ~isnumeric(digital)
    if strcmp(digital,'fs')
        digmax = 1;
    else
        error('sorry, enter either "fs" or a bit rate');
    end
elseif isnumeric(digital)
        digmax = 2^(digital - 1);
end

   
gains = 1:1:11;
vpk2pk = [1.61,1.573351926,1.370332242,...
1.05154019,0.678931371,0.509126703,...
0.344211896,0.212239335,0.130865713,...
0.056471052,0.021969788];
vmax = vpk2pk./2;
sensitivity = 715.7574625; 

if gain > 0
idx = find(gains == gain);
vpk = [-vmax(idx),vmax(idx)];
end

if strcmp(analog,'mV') && gain > 0
    factor = (vmax(idx)*1000)/digmax;
elseif strcmp(analog,'V') && gain > 0
    factor = vmax(idx)/digmax;
elseif strcmp(analog,'Pa') && gain > 0
    factor = ((vmax(idx)*1000)/sensitivity)/digmax;
elseif gain == 0
    factor = 1.0;
    vpk = [-digmax,digmax-1];
else
    error('Sorry what? Enter either mV, V, or Pa');
end

end