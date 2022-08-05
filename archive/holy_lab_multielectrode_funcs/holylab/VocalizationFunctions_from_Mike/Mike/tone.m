function y = tone(freq, t, fs, maxamp)

rfreq = freq*2*pi;
times = linspace(0,t,fs*t);

y = maxamp*cos(rfreq*times);


end
