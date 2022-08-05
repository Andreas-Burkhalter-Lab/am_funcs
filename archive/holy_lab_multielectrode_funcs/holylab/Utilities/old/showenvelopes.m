envfiles = dirbyname('*.env');
ephysbase = ephysfromai(envfiles{1});
ephysbase.envelopefile = envfiles{1};
ephysbase.tag = '';
ephys5 = ephyssubrange(ephysbase,[1 50000]);
ephys5 = ephysfetch(ephys5,'envelope');
epp = struct('tags','','fieldtoplot','envelope');
for i = 1:64;
  subplot(8,8,i)
  epp.number = i;
  ephysplot(ephys5,epp)
  disp(num2str(ephysbase.channels(i)))
  title(num2str(ephysbase.channels(i)));
end