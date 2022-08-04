function [ listout ] = RD_lists_traces2average( whichlist )
%RD_LISTS_TRACES2AVERAGE Lists of .xsg trace files to be combined into an
%average. 
%   Input argument 'whichlist' should be a string matching a list name below. 
%   Lists are two-column cell arrays, with the first column specifying the file
%   and the second column specifying the trace (1 or 2) to be analyzed. 
%   Each list should be called by a separate switch-case specified 
%   by the list name.  This function is called by RD_avg_xsg.
%%%% Last updated 6/15/15 flash drive

%% Lists
%%% First column = filename, second column = trace number within this file.
switch whichlist
    case 'samplelist'
        listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\110615\RD0009\RD0009AAAA0015.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\110615\RD0009\RD0009AAAA0016.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\110615\RD0009\RD0009AAAA0017.xsg', 1;
        };
    
    case 'pair_150615_RD0010_100pA'
        listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0002.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0007.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0012.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0017.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0022.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0027.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0032.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0037.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0042.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0047.xsg', 1;
        };
        
    case 'pair_150615_RD0010_200pA'
        listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0003.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0008.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0013.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0018.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0023.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0028.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0033.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0038.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0043.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0048.xsg', 1;
        };
    
    case 'pair_150615_RD0010_300pA'
        listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0004.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0009.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0014.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0019.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0024.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0029.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0034.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0039.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0044.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0049.xsg', 1;
        };
    
    case 'pair_150615_RD0010_400pA'
          listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0005.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0010.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0015.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0020.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0025.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0030.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0035.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0040.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0045.xsg', 1;
        };
    
     case 'pair_240615_RD0003_300pA'
         listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0004.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0009.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0014.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0019.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0024.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0029.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0034.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0039.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0044.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0049.xsg', 2;
        };
    
    
    case 'pair_240615_RD0003_400pA'
         listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0005.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0010.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0015.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0020.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0025.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0030.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0035.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0040.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0045.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0050.xsg', 2;
        };
    
    case 'L23_allpairs_100pA'
        listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0002.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0007.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0012.xsg', 1;
%        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0017.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0022.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0027.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0032.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0037.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0042.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0047.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0002.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0007.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0012.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0017.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0022.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0027.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0032.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0037.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0042.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0047.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0002.xsg', 1;
%        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0007.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0012.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0017.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0022.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0027.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0032.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0037.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0042.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0047.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0017\RD0017AAAA0002.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0017\RD0017AAAA0007.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0017\RD0017AAAA0012.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0017\RD0017AAAA0017.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0018\RD0018AAAA0002.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0018\RD0018AAAA0007.xsg', 1;
 %       'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\130715\RD0004\RD0004AAAA0002.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\130715\RD0004\RD0004AAAA0007.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0002.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0007.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0012.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0022.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0027.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0002.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0012.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0022.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0027.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0032.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0002.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0007.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0012.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0022.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0027.xsg', 1;
        };
        
     case 'L23_allpairs_200pA'    
        listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0003.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0008.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0013.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0018.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0023.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0028.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0033.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0038.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0043.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0048.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0003.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0008.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0013.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0018.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0023.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0028.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0033.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0038.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0043.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0048.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0003.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0008.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0013.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0018.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0023.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0028.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0033.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0038.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0043.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0048.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0017\RD0017AAAA0003.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0017\RD0017AAAA0008.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0017\RD0017AAAA0013.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0017\RD0017AAAA0018.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0018\RD0018AAAA0003.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0018\RD0018AAAA0008.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\130715\RD0004\RD0004AAAA0003.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\130715\RD0004\RD0004AAAA0008.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0003.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0008.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0013.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0023.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0028.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0003.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0018.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0023.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0028.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0033.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0003.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0008.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0013.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0023.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0028.xsg', 1;
        };
    
      case 'L23_allpairs_300pA'    
        listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0004.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0009.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0014.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0019.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0024.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0029.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0034.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0039.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0044.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0049.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0004.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0009.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0014.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0019.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0024.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0029.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0034.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0039.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0044.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0049.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0004.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0009.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0014.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0019.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0024.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0029.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0034.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0039.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0044.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0049.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0017\RD0017AAAA0004.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0017\RD0017AAAA0009.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0017\RD0017AAAA0014.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0017\RD0017AAAA0019.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0018\RD0018AAAA0004.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0018\RD0018AAAA0009.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\130715\RD0004\RD0004AAAA0004.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\130715\RD0004\RD0004AAAA0009.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0004.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0009.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0014.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0024.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0029.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0004.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0019.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0024.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0029.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0034.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0004.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0009.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0014.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0024.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0029.xsg', 1;
        };
    
      case 'L23_allpairs_400pA'    
        listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0005.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0010.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0015.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0020.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0025.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0030.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0035.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0040.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0045.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\150615\RD0010\RD0010AAAA0050.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0005.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0010.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0015.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0020.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0025.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0030.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0035.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0040.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0045.xsg', 2;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\240615\RD0003\RD0003AAAA0050.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0005.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0010.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0015.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0020.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0025.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0030.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0035.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0040.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0044.xsg', 1;
%         'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0006\RD0006AAAA0050.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0017\RD0017AAAA0005.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0017\RD0017AAAA0010.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0017\RD0017AAAA0015.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0017\RD0017AAAA0020.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\June\250615\RD0018\RD0018AAAA0010.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\130715\RD0004\RD0004AAAA0005.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\130715\RD0004\RD0004AAAA0010.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0005.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0010.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0015.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0025.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0003\RD0003AAAA0030.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0005.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0020.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0025.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0030.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0008\RD0008AAAA0035.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0005.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0010.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0015.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0025.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\140715\RD0013\RD0013AAAA0030.xsg', 1;
        };
    
    case 'L5_allpairs_100pA'
        listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0003.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0008.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0013.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0023.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0033.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0002.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0007.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0012.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0017.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0022.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0002.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0007.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0012.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0017.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0022.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0002.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0007.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0012.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0017.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0022.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0002.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0007.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0012.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0022.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0002.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0007.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0012.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0017.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0022.xsg', 2;
        };
    
    case 'L5_allpairs_200pA'
        listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0004.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0009.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0014.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0029.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0034.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0003.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0008.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0013.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0018.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0023.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0003.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0008.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0013.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0018.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0023.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0003.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0008.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0013.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0018.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0023.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0003.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0008.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0013.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0018.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0023.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0003.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0008.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0013.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0018.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0023.xsg', 2;
            };
      
    case 'L5_allpairs_300pA'
        listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0005.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0010.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0015.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0030.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0035.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0004.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0009.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0014.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0019.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0024.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0004.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0009.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0014.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0019.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0024.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0004.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0009.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0014.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0019.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0024.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0004.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0009.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0014.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0019.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0024.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0004.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0009.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0014.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0019.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0024.xsg', 2;
            };
    
     case 'L5_allpairs_400pA'
        listout = {
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0006.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0011.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0016.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0031.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\July\240715\RD0004\RD0004AAAA0036.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0005.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0010.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0015.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0020.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\140815\RD0002\RD0002AAAA0025.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0005.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0010.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0015.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0020.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0003\RD0003AAAA0025.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0005.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0010.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0015.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0020.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\August\170815\RD0011\RD0011AAAA0025.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0005.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0010.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0015.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0020.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0002\RD0002AAAA0025.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0005.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0010.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0015.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0020.xsg', 2;
        'C:\Ephus\R210\Ephus\Data\RinaldoData\2015\September\030915\RD0009\RD0009AAAA0025.xsg', 2;
            };
        
    case 'sampleg'
        listout = {
        'G:\RD0009\RD0009AAAA0015.xsg', 1;
        'G:\RD0009\RD0009AAAA0016.xsg', 1;
        'G:\RD0009\RD0009AAAA0017.xsg', 1;
        };
end



end

