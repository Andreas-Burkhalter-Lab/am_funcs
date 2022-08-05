gIsRecording = 1;
hax = gca;
htitle = title('Starting...');
nfreq = 256;
navg = 2;
scanrate = 200000;
npix = 256;
nperblock = 20;
nblocks = 5*16;
shg
sng = RecordingWhis(2,scanrate,nfreq,navg,nperblock,nblocks,1,npix,hax,htitle,'Test save','save.bin');
