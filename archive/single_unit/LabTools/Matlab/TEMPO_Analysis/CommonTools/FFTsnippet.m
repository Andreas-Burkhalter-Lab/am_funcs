FT = abs(fft(g - mean(g), 128))
FT2 = FT(1:length(FT)/2)
[max_ampl max_indx] = max(FT2);
plot(FT2)