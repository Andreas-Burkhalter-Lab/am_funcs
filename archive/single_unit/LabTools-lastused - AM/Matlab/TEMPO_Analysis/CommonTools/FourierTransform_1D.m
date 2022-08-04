%-----------------------------------------------------------------------------------------
%-- FourierTransform_1D.m: This function takes in an (x,y) pair of data vectors (data_x, data_y)
%--     and computes the Fourier amplitude spectrum.  FFT_PTS is the number of points used in
%--     computing the DFT; DC_remove should be set to 1 if you want the DC component removed, and
%--     plot_flag should be set to 1 if you want the results to be plotted for you.
function [freq, ampl] = FourierTransform_1D(data_x, data_y, FFT_PTS, DC_remove, plot_flag)

DFT_input = data_y - (DC_remove == 1)*mean(data_y);
DFT_mag = abs(fft(DFT_input, FFT_PTS)); %abs computes the complex modulus when input is complex
half_length = floor(length(DFT_mag)/2);
DFT_mag_half = DFT_mag(1:half_length);
ampl = DFT_mag_half;

dx = (data_x(length(data_x)) - data_x(1))/(length(data_x)-1);
max_freq = 1/(2*dx);
freq = (0:(length(DFT_mag_half)-1)).*(max_freq/length(DFT_mag_half));

if (plot_flag)
    figure;
    subplot(2,1,1);
    bar(data_x, data_y, 'k-');
    title('Input Data');
    
    subplot(2,1,2);
    plot(freq, DFT_mag_half, 'k-');
    title('Amplitude Spectrum');
end

return;
