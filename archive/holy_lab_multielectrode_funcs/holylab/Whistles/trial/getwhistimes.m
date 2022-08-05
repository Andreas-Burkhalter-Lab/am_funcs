bp = struct('slope',1.25,'widthp',8,'widthm',8,'aratio',1,'npix',129,'freqrange',[60 110])
dbar = buildbar(bp);
convo = sngconv(sng,dbar);
tw = whistimes(convo,struct('thresh',20,'width',150));
%clf
%hax(1) = subplot(3,1,1);
%plot(convo);
%hax(2) = subplot(3,1,2);
%plot([tw;tw],[-1 1],'r')
%set(gca,'XLim',[1 size(sng,2)])
%hax(3) = subplot(3,1,3);
%spsngplot(sng)
%SliderWindow(hax([3 1 2]),struct('axisupdatefcn',[@setspimg @set @set]))
['Found ',num2str(length(tw)),' whistles']
snips = whissnipmem(sng,tw-30,150,[1 129]);
return
snips16 = whissnipmem(sng,tw-30,150,[1 129],[16 16]);
figure
for i = 1:size(snips,3); imagesc(abs(snips(:,:,i))); axis xy; title(num2str(i)); pause; end