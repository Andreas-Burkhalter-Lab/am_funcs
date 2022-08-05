%[waveform] = createwaveform(1000,{[200],[400],[700]},{[400 980 616 -132],[500 1000 500 -200],[1000 500]});
snipform = [1 4 9 16 18 19.5 20 19.5 18 16 9 4 1]*1000;
[waveform] = createwaveform(1100,[100 200.1 300.2 400.3 500.4 600.5 700.6 800.7 900.8 1000.9]-(length(snipform)+1)/2,snipform);
scanrange = [1 length(waveform)-1];
channels = [1];
thresh = [1 ];
sniprange = [-10 30];


[tsnip,snip] = snippetfile(waveform,scanrange,channels,thresh,sniprange,struct('interptimes',1,'interpsnips',1));
