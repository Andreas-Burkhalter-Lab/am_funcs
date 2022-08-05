for i = 0:63
  [snip,time,h] = LoadSnip('01spont.ssnp',i);
  if ~isempty(snip)
    plot(snip)
    title(['Channel ' num2str(i)])
    shg
    pause
  end
end
