function bar = buildbar(barparms)
%   barparms is a structure specifying the bar, with fields
%     slope: determines the tilt of the bar, in freq/time (pixel units, not physical units).  This should be matched to the slope of the whistle.
%     widthp: the width, in pixels, of the center of the bar along the frequency axis.
%     widthm: the width, in pixels, of the negative-going side lobes along the frequency axis.
%     aratio: the absolute value of the ratio between the amplitude of the lobe to the center
%     npix: the total number of pixels along the frequency axis
%     freqrange: the frequency range (in pixels) of support for the bar center
bwidth = ceil((diff(barparms.freqrange)-barparms.widthp)/barparms.slope);
bar = zeros(barparms.npix,bwidth);
for i = 0:bwidth-1
bot = floor(i*barparms.slope + barparms.freqrange(1));
bar(bot:bot+barparms.widthp,i+1) = 1;
bar(max(1,bot-barparms.widthm):bot-1,i+1) = -barparms.aratio;
bar(min(barparms.npix,bot+barparms.widthp+1):min(barparms.npix,bot+barparms.widthp+barparms.widthm),i+1) = -barparms.aratio;
end
bar = sparse(bar);

