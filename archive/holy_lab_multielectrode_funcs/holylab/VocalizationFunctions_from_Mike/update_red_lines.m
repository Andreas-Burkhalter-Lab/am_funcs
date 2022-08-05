function update_red_lines(twhis)

newXaxis= get(gca,'xlim')
these_ones=find(twhis(1,:)>=newXaxis(1)& twhis(2,:)<=newXaxis(2));
for x=these_ones
  plot(twhis(:,x),[30 30],'r')
end