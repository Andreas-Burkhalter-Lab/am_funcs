function listout = RD_lists_spikeLatency( whichlist)
%RD_LISTS_SPIKELATENCY Lists of .xsg trace files to analyze for
%relationship between presynaptic input vs. postsynaptic first-spike
%latency. 

% The number following each file is the channel (1 or 2) which will be considered
% the postsynaptic (pyramidal) cell; we will look in this channel to for
% the timing of the first spike. The onset of stimulus to the other channel 
% (1 or 2) will be used to determine the latency of first spike to stimulus
% onset.

%% Lists
%%% First column = filename, second column = trace number of postsynaptic neuron within this file.
switch whichlist
    case 'samplelist'
        listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0004.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0005.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0006.xsg', 2;
        };
    
    case 'L5_allpairs_0pA'
        listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0001.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0006.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0011.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0016.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0021.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0026.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0031.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0036.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0041.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0046.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0002.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0007.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0012.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0017.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0022.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0027.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0032.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0037.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0042.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0047.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0001.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0006.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0011.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0016.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0021.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0026.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0031.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0036.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0041.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0046.xsg', 2;
        };
      
    case 'L5_allpairs_100pA'
        listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0002.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0007.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0012.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0017.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0022.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0027.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0032.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0037.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0042.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0047.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0003.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0008.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0013.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0018.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0023.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0028.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0033.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0038.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0043.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0048.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0002.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0007.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0012.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0017.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0022.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0027.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0032.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0037.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0042.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0047.xsg', 2;
        };
    
    case 'L5_allpairs_200pA'
        listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0003.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0008.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0013.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0018.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0023.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0028.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0033.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0038.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0043.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0048.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0004.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0009.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0014.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0019.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0024.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0029.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0034.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0039.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0044.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0049.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0003.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0008.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0013.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0018.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0023.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0028.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0033.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0038.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0043.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0048.xsg', 2;
        };
      
    case 'L5_allpairs_300pA'
        listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0004.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0009.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0014.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0019.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0024.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0029.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0034.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0039.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0044.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0049.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0005.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0010.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0015.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0020.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0025.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0030.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0035.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0040.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0045.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0050.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0004.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0009.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0014.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0019.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0024.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0029.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0034.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0039.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0044.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0049.xsg', 2;
        };
    
     case 'L5_allpairs_400pA'
        listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0005.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0010.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0015.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0020.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0025.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0030.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0035.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0040.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0045.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0009\RD0009AAAA0050.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0006.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0011.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0016.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0021.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0026.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0031.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0036.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0041.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0046.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0004\RD0004AAAA0001.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0005.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0010.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0015.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0020.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0025.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0030.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0035.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0040.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0045.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0008\RD0008AAAA0050.xsg', 2;
            };
        
    case 'sampleg'
        listout = {
        'G:\RD0009\RD0009AAAA0015.xsg', 1;
        'G:\RD0009\RD0009AAAA0016.xsg', 1;
        'G:\RD0009\RD0009AAAA0017.xsg', 1;
        };
end



end