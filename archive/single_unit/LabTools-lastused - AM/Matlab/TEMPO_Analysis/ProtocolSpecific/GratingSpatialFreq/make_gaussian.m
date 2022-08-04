%This program was used to make a gaussian function to help eliminate high
%frequency noise from fourier analyses

%GAUSSFUNC Used by GAUSSFIT.
%  GAUSSFUNC assumes a function of the form
%
%	  y = q(1) + q(2) * exp(-0.5*((x - q(3))/q(4)).^2 )
%
%	thus q(1) is base rate, q(2) is amplitude, q(3) is center, and q(4) is size
sigma = 0.1;
dt = 0.005;
x = [-0.3:dt:0.3];
q = [0, 1, 0, sigma];
z = q(1) + q(2) * exp(-0.5*((x - q(3))/ q(4)).^2);
plot(x,z)



t = 0:0.001:0.6;
x = sin(2*pi*50*t)+sin(2*pi*120*t);
y = x + 2*randn(size(t));
plot(1000*t(1:50),y(1:50))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('time (milliseconds)')

Y = fft(y,512)
Pyy = Y.* conj(Y) / 512;
f = 1000*(0:256)/512;
plot(f,Pyy(1:257))
title('Frequency content of y')
xlabel('frequency (Hz)')