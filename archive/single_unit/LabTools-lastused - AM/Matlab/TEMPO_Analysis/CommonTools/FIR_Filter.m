% FIR_Filter.m: This function implements a finite impulse response high/low pass filter.
% Greg DeAngelis, 4/22/03
% x is the input vector to be filtered
% CutoffFreq is the desired cutoff frequency in Hz
% SampleFreq is the sampling frequency of the input data in Hz
% Type determines the filter type: 'low' for LowPass, 'high' for HighPass
% Order is the order of the filter -- use a higher integer for a steeper rolloff (must be >2)
% BodePlot is a flag (0 or 1) that determines whether the frequency response plot of the filter will
%   be graphed

function [y] = FIR_Filter(x, CutoffFreq, SampleFreq, Type, Order, BodePlot)

if (CutoffFreq > SampleFreq/2)
    disp('Error: The cutoff frequency must be less than half the sampling frequency');
    return;
end

frac = CutoffFreq/(SampleFreq/2);

[Bb, Aa] = fir1(Order, frac, Type);

if (BodePlot)
    freqz(Bb, Aa);
end

y = filter(Bb, Aa, x);

if (BodePlot)
    t = 0:(length(x)-1);
    t = t * (1/SampleFreq);
    [ampl, freq] = FourierTransform_1D(t, y, length(y), 0, 1);
end

return;