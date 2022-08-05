function [snip,f,t] = whissnipmem(sng,tsnip,twidth,band,sizeout)
if (nargin < 5)
        sizeout = [];
end
if (nargin < 4 | isempty(band))
        band = [35 115];
end
nsnip = length(tsnip);
nfft = 256;
f = 1:size(sng,1);
fsi = find(f >= band(1) & f < band(2));
nfreq = length(fsi);
%nt = floor(2*twidth*h.scanrate/nfft)-1;
nt = twidth;
if (~isempty(sizeout))
        tbin = linspace(1,nt,sizeout(1)+1);
        tbinf = floor(tbin);
        tbinc = ceil(tbin);
        fbin = linspace(1,nfreq,sizeout(2)+1);
        fbinf = floor(fbin);
        fbinc = ceil(fbin);
        snip = zeros(sizeout(1),sizeout(2),nsnip);
else
        snip = zeros(nfreq,nt,nsnip);
end
for k = 1:nsnip
        %z = loadmc(filename,tsnip(k)+[0 twidth]);
        %B = specgram(z,nfft);
    B = sng(:,tsnip(k):tsnip(k)+twidth-1);
        % Now do the decimation
        if (~isempty(sizeout))
                temp = abs(B(fsi,:));
                for i = 1:sizeout(1)
                        for j = 1:sizeout(2)
                                % Sum over the whole block
                                tmpsum = sum(sum(temp(fbinc(i):fbinf(i+1),tbinc(j):tbinf(j+1))));
                                % Now pick up the fractional parts on the edges
                                tmpsum = tmpsum + (fbinc(i)-fbin(i))*sum(temp(fbinf(i),tbinc(j):tbinf(j+1)));
                                tmpsum = tmpsum + (fbin(i+1)-fbinf(i+1))*sum(temp(fbinc(i+1),tbinc(j):tbinf(j+1)));
                                tmpsum = tmpsum + (tbinc(j)-tbin(j))*sum(temp(fbinc(i):fbinf(i+1),tbinf(j)));
                                tmpsum = tmpsum + (tbin(j+1)-tbinf(j+1))*sum(temp(fbinc(i):fbinf(i+1),tbinc(j+1)));
                                % Now do the corners
                                tmpsum = tmpsum + (fbinc(i)-fbin(i))*(tbinc(j)-tbin(j))*temp(fbinf(i),tbinf(j));
                                tmpsum = tmpsum + (fbinc(i)-fbin(i))*(tbin(j+1)-tbinf(j+1))*temp(fbinf(i),tbinc(j+1));
                                tmpsum = tmpsum + (fbin(i+1)-fbinf(i+1))*(tbinc(j)-tbin(j))*temp(fbinc(i+1),tbinf(j));
                                tmpsum = tmpsum + (fbin(i+1)-fbinf(i+1))*(tbinc(j)-tbin(j))*temp(fbinc(i+1),tbinc(j+1));
                                % Whew!
                                snip(i,j,k) = tmpsum;
                        end
                end
        else
                snip(:,:,k) = B(fsi,:);
        end
end
f = f(fsi);
t = linspace(0,twidth,nt);
