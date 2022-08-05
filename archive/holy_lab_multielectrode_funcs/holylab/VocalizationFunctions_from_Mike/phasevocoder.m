function y = phasevocoder(varargin)
% y = phasevocoder(x,nfft,divfac,thresh,bandsig)
% y = phasevocoder(filename,trange,nfft,divfac,thresh,bandsig)
% bandsig is measured in normalized units (relative to fc), not
% physical units
% Copyright 2001 by Timothy E. Holy <holy@pcg.wustl.edu>
error(nargchk(3,6,nargin))
[nfft,divfac,thresh,npts,fromfile,x,fid,bandsig] = pvargchk(varargin);
nffts = nfft/divfac;
% Fourier window stuff
wa = hanning(nfft)';
ws = hanning(nfft/divfac)';
noverlapa = ceil(3*nfft/4);
noverlaps = ceil(3*nfft/4/divfac);
%noverlapa = ceil(7*nfft/8);
%noverlaps = ceil(7*nfft/8/divfac);
wsn = windownorm(ws,noverlaps);
% Compute sizes, number of blocks, etc
skipa = nfft-noverlapa;
skips = nffts - noverlaps;
%[nfft nffts skipa skips]
nblocks = floor((npts-nfft)/skipa) + 1;
nptsout = skips*nblocks+nffts;
y = zeros(1,nptsout);
ts = 0:nffts-1;
fa = 2*pi*(0:nfft-1)/nfft;
fs = 2*pi*(0:nffts-1)/nffts;
if (~isempty(bandsig))
        findx = find(fa > pi*bandsig(1) & fa < pi*bandsig(2));
else
        findx = 1:nfft;
end
%findx
%angckeep = zeros(nblocks,nfft);
%B = specgram(x,nfft,[],nfft,noverlapa);
%angB = angle(B);
%angcorr = unwrap(angB');
wvout = zeros(0,nffts);
%nblocks = 5;
%npts
for k = 1:nblocks
        if (fromfile)
                z = ReadBinaryData(fid,1,[0 nfft-1]);
                % Now position for the next read
                if (k < nblocks)
                        status = fseek(fid,-noverlapa*2,'cof');        % *2 bcs int16
                        if (status)
                                error(ferror(fid))
                        end
                end
        else
                z = x((k-1)*skipa + (1:nfft));
        end
        % Compute windowed transform, then threshold
        zt = fft(z.*wa);
        abszt = abs(zt);
        %semilogy(abszt);
        %semilogy([1 length(abszt)],[thresh thresh],'k--');
        %drawnow
        ikeep = find(abszt(findx) > thresh);
        ikeep = findx(ikeep);
        % Compute the phase derivative
        ang = angle(zt);
        if (k > 1)
                dang = diff([anglast;ang]) - fa*skipa;        % Phase difference - expected drift
                dang = mod(dang+pi,2*pi)-pi;
                %dang = mod(dang,2*pi);
                dang([1 (nfft/2+1) end]) = 0;        % The real components have phase 0
        else
                angclast = ang;
                dang = zeros(1,nfft);
        end
        anglast = ang;
        % Now synthesize the frequency-divided waves
        winout = zeros(1,nffts);
        omega = fa + dang/skipa;        % This is the true instantaneous frequency
        %angclast(2:2:end) = angclast(2:2:end)+pi;
        %omega = fa;
        %angclast = ang;
        % Control output with the following:
        if (k < 1)
                fi = 24:29;
                'Start'
                dang(fi)
                fa(fi)
                'Omega'
                omega(fi)
                angclast(fi)
                ang(fi)
                mod(angclast(fi)-ang(fi)+pi,2*pi)-pi
        end
        %omega = fs;
        %ang0 = ang;
        for i = 1:length(ikeep)
                %winout = winout + abszt(ikeep(i)) * cos( fs(ikeep(i))*ts + dangc(ikeep(i))/(divfac*skipa) );
                %winout = winout + abszt(ikeep(i)) * cos( (fs(ikeep(i)) + dangc(ikeep(i))/(divfac*skipa))*ts );
                wintmp = abszt(ikeep(i)) * cos(omega(ikeep(i))*ts + angclast(ikeep(i)));
                winout = winout + wintmp;
        end
        winout = winout ./ wsn;        % Correct the normalization
        % Add the synthesized waveform to the total accumulated waveform
        indx = (k-1)*skips + (1:nffts);
        y(indx) = y(indx) + winout/nffts;
        % The next lines seem strange, comment them out (TEH 2004-08-25)
        %if (k < 5)
        %        wvout(end+1,:) = abszt(fi(1)) * cos(omega(fi(1))*ts + angclast(fi(1)));
        %end
        angclast = angclast+omega*skips;                        % Update the phase offset
end
%keyboard


function [nfft,divfac,thresh,npts,fromfile,x,fid,bandsig] = pvargchk(P);
% y = phasevocoder(x,nfft,divfac,thresh)
% y = phasevocoder(filename,trange,nfft,divfac,thresh)
if (ischar(P{1}))
        filename = P{1};
        trange = P{2};
        nfft = P{3};
        divfac = P{4};
        if (length(P) >= 5 & ~isempty(P{5}))
                thresh = P{5};
        else
                thresh = 0;
        end
        if (length(P) >= 6)
                bandsig = P{6};
        else
                bandsig = [];
        end
        fromfile = 1;
        
        [fid,npts] = posaifile(filename,trange);
        x = [];
else
        x = P{1};
        nfft = P{2};
        divfac = P{3};
        if (length(P) >= 4 & ~isempty(P{4}))
                thresh = P{4};
        else
                thresh = 0;
        end
        if (length(P) >= 5)
                bandsig = P{5};
        else
                bandsig = [];
        end
        npts = length(x);
        fromfile = 0;
        fid = [];
end
nffto = nfft;
nfft = 2^ceil(log2(nfft));
if (nfft ~= nffto)
        warning(sprintf('nfft changed to %d',nfft));
end
divfaco = divfac;
divfac = 2^floor(log2(divfac));
if (divfac ~= divfaco)
        warning(sprintf('divfac changed to %d',divfac));
end
if (nfft <  4)
        error('nfft must be >= 4!');
end
