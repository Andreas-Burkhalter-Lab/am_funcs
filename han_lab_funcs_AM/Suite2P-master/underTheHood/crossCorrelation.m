function [ymax, xmax, cMax] = crossCorrelation(T, I, phasecorr)

I = single(I);
T = single(T)/100;

% compute average power spectrum and divide by that
T = my_conv2(T, 1, [1 2]);
f2T = conj(fft2(T));

mNorm = abs(f2T);

% f2T = f2T./mNorm;

ccF = f2T .* fft(fft(I, [],1),[],2);

if nargin>2 && ~isempty(phasecorr)
    ccF = ccF./(1e-2 + abs(ccF));
end
cc  = ifft(ifft(ccF, [], 1), [], 2);

cc = fftshift(fftshift(cc, 1), 2);
cc = real(cc);

cc  = my_conv2(cc, 0.5, [1 2]);


[cMax, iMax] = max(reshape(cc, [], size(cc,3)), [],1);

[Ly, Lx, NT] = size(cc);
ymax = rem(iMax-1, Ly) + 1;
xmax = ceil(iMax/Ly);


ymax = ymax - ceil((Ly+1)/2);
xmax = xmax - ceil((Lx+1)/2);