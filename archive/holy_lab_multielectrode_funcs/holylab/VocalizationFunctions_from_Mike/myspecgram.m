function [B,ang,angc] = myspecgram(x,nfft)
npts = length(x);
noverlap = ceil(0.75*nfft);
skip = nfft-noverlap;
nblocks = floor((npts-nfft)/skip) + 1;
wa = hanning(nfft)';
zt = zeros(1,nfft);
B = zeros(nblocks,nfft);
ang = zeros(nblocks,nfft);
angc = zeros(nblocks,nfft);
angclast = zt;
for k = 1:nblocks
        z = x((k-1)*skip + (1:nfft));
        zt = fft(z.*wa);
        B(k,:) = zt;
        angtemp = angle(zt);
        ang(k,:) = angtemp;
        dang = diff([angclast;angtemp]);
        dang = mod(dang+pi,2*pi)-pi;
        angnew = angclast + dang;
        if (k == 2)
                %size(angtemp)
                %size(angctemp)
                %angtemp(:,1:8)
                [angclast(1:8);angtemp(1:8)]
                [angclast(1:8);angnew(1:8)]
                %angclast(1:8)
        end
        angc(k,:) = angnew;
        angclast = angnew;
end
