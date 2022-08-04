%LowPassFilter_DiffEqn.m: This function implements a low-pass filter (one pole)
%  as a difference equation.  The input, x, is a time sequence of numbers, and y is the filtered output.  
%  CutFreq is the cutoff frequency (-3dB point) (in Hz).
%  Fsamp is the sampling frequency of the input data (in Hz)

%NOTE: if you want a sharper filter, you can apply this recursively multiple times to your data.
function [y] = LowPassFilter_DiffEqn(x, CutFreq, Fsamp);

Pole =1/(2*pi*CutFreq);

PoleTerm = Pole*2*Fsamp;

y(1) = 0;
for i=2:length(x)
    y(i) = (1/(1+PoleTerm))*(-y(i-1)*(1-PoleTerm) + x(i) + x(i-1));
end

return;