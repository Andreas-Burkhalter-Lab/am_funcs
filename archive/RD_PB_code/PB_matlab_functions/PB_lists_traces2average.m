function [ listout ] = PB_lists_traces2average( whichlist )
%PB_LISTS_TRACES2AVERAGE Lists of .xsg trace files to be combined into an
%average. 
%   Input argument 'whichlist' should be a string matching a list name below. 
%   Lists are two-column cell arrays, with the first column specifying the file
%   and the second column specifying the trace (1 or 2) to be analyzed. 
%   Each list should be called by a separate switch-case specified 
%   by the list name.  This function is called by RD_avg_xsg.
%%%% Last updated 6/19/15 ephus comp

%% Lists
%%% First column = filename, second column = trace number within this file.
switch whichlist
    
     case '050515_PB0002_Patch'
       listout = {
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0074.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0075.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0077.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0078.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0079.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0080.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0084.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0085.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0086.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0087.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0088.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0089.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0091.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0092.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0093.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0094.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0095.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0097.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0098.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0100.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0103.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0107.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0109.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0110.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0113.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0115.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0117.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0118.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0121.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0123.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0125.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0129.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0131.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0132.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0139.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0141.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0143.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0144.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0145.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0146.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0149.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050515\PB0002\Paired recording\PB0002\PB0002AAAA0150.xsg', 1;
       };
    
   case '050615_PB0001_Patch'
       listout = {
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0293.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0294.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0295.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0298.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0299.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0300.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0303.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0304.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0306.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0309.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0310.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0312.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0314.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0316.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0318.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0324.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0327.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0328.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0330.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0333.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0334.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0335.xsg', 1;
       'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\050615\PB0001\Paired recording\PV2Pyr\PB0001\Trial 2\PB0001\PB0001AAAA0335.xsg', 1;
       };
    
 
     
     
    case '052815_PB0001_Patch'
       listout = {
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0091.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0096.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0099.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0104.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0110.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0114.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0115.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0119.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0124.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0130.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0132.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0135.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0137.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0140.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0142.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0143.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0145.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0148.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0150.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0151.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0154.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0156.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0162.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0165.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0172.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0176.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0178.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0180.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0187.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0193.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0194.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0204.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0205.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0212.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0219.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0226.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0227.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0001\After Blocker\PB0001AAAA0228.xsg', 1;
        
        };
        
    case '052815_PB0002_Patch'
       listout = {
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0025.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0026.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0054.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0055.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0063.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0065.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0066.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0067.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0068.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0069.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0070.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0071.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0072.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0073.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0074.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0075.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0076.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0077.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0078.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0079.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0080.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0081.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0082.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0083.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0084.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0085.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0087.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0088.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0089.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0090.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0091.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0092.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0093.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0094.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0095.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0096.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0097.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0098.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0100.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0101.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0102.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0103.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0104.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0105.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052815\PB0002\PV2Pyr\PB0002AAAA0108.xsg', 1;
         };
            
    case '052915_PB0006_Patch'
        listout = {
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0031.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0032.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0033.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0034.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0035.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0036.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0037.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0038.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0039.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0040.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0041.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0042.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0043.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0044.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0045.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0046.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0047.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0048.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0049.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0050.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0051.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0052.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0053.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0054.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0055.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0056.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0057.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0058.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0059.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0060.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0061.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0062.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0063.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0064.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0065.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0066.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0067.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0068.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0069.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0070.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0072.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0073.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0074.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0075.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0076.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0077.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0078.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0079.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0080.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0084.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0089.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0091.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0095.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0096.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0098.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0006\PB0006AAAA0099.xsg', 1;
         };
        
    case '052915_PB0007_Interpatch'
        listout = {
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0016.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0017.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0018.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0020.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0021.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0024.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0025.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0026.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0027.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0029.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0034.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0036.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0037.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0038.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0039.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0040.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0042.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0045.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0047.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0048.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0050.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0051.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0052.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0053.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0054.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0055.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0057.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0058.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0059.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0061.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0062.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0063.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0064.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0065.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0066.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0067.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0068.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0069.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0072.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0074.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0075.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0076.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0078.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0080.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0082.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0083.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0084.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0085.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0086.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0087.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0088.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0089.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0090.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0091.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0092.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0093.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0094.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0095.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\052915\PB0007\PB0007AAAA0097.xsg', 1;
         };
    
    
    
    case '061215_PB0003_Interpatch'
        listout = {
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0012.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0013.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0014.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0015.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0016.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0017.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0018.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0019.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0020.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0021.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0024.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0026.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0028.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0032.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0035.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0036.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0038.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0039.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0040.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0041.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0042.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0043.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0045.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0046.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0049.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0051.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0057.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0003\PB0003AAAA0058.xsg', 1;
        
        };
   
    case '061215_PB0004_patch'
        listout = {
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0066.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0069.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0070.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0071.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0077.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0078.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0079.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0080.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0081.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0082.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0088.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0089.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0090.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0093.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0096.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0097.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0102.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0103.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0104.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0105.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0110.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0111.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0112.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0114.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0117.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0118.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0122.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0132.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0135.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0140.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061215\PB0004\PB0004AAAA0155.xsg', 1;
       
        };   
        
   
    case '061915_PB0006_Interpatch'
        listout = {
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0029.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0031.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0032.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0036.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0037.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0038.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0065.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0091.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0092.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0093.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0099.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0100.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0103.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0104.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0113.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0114.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0119.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0120.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0156.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0161.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0168.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0170.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0171.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0174.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0175.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0176.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0180.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0181.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0182.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0183.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0184.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0186.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0187.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0189.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0190.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0194.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0199.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0206.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0006\PB0006AAAA0211.xsg', 1; 
        };
    
    
    case '061915_PB0005_patch'
         listout = {
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0005\PB0005AAAA0075.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0005\PB0005AAAA0089.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0005\PB0005AAAA0098.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0005\PB0005AAAA0118.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0005\PB0005AAAA0131.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0005\PB0005AAAA0187.xsg', 1;
           };
       
    case '061915_PB0001_patch'
         listout = {
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0001\PB0001AAAA0164.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0001\PB0001AAAA0192.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0001\PB0001AAAA0230.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0001\PB0001AAAA0247.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\061915\PB0001\PB0001AAAA0265.xsg', 1;
           };
       
    case '062015_PB0001_Interpatch'
         listout = {
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0134.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0135.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0137.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0138.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0139.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0140.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0141.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0142.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0143.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0144.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0145.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0146.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0147.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0149.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0150.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0151.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0152.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0153.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0154.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0155.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0156.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0157.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0158.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0159.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0160.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0162.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0163.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0166.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0167.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0172.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0177.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0178.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0180.xsg', 1
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0181.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0185.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0186.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0188.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0189.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062015\PB0001\With DNQX and CPP\PB0001AAAA0190.xsg', 1;
        
           };
       
     case '062315_PB0003_Patch'
         listout = {  
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062315\PB0005\PB0005AAAA0001.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062315\PB0005\PB0005AAAA0002.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062315\PB0005\PB0005AAAA0004.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062315\PB0005\PB0005AAAA0005.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062315\PB0005\PB0005AAAA0008.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062315\PB0005\PB0005AAAA0011.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062315\PB0005\PB0005AAAA0012.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062315\PB0005\PB0005AAAA0017.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062315\PB0005\PB0005AAAA0046.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062315\PB0005\PB0005AAAA0049.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062315\PB0005\PB0005AAAA0060.xsg', 1;
        
         };
     
     case '062615_PB0003_Patch_Before Picrotoxin'
         listout = { 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0011.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0012.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0014.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0015.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0016.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0017.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0018.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0019.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0020.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0021.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0022.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0023.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0024.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0025.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0026.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0027.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0028.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0029.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0030.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0031.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0032.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0033.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0034.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0035.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0036.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0037.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0038.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0039.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0040.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0041.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0043.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0044.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0045.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0046.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0047.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0048.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\Before Picrotoxin\PB0003AAAA0049.xsg', 1;
        };
        
      case '062615_PB0003_Patch_After Picrotoxin'
         listout = { 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\After Picrotoxin\PB0003AAAA0091.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\After Picrotoxin\PB0003AAAA0092.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\After Picrotoxin\PB0003AAAA0093.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\After Picrotoxin\PB0003AAAA0094.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\After Picrotoxin\PB0003AAAA0095.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\After Picrotoxin\PB0003AAAA0096.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\After Picrotoxin\PB0003AAAA0097.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\After Picrotoxin\PB0003AAAA0098.xsg', 1;
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062615\PB0003\After Picrotoxin\PB0003AAAA0099.xsg', 1;
        
     };
      
     case '062915_PB0002_Interpatch_Before Picrotoxin_IPSC'
         listout = { 
             
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0001.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0002.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0003.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0004.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0005.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0006.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0007.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0008.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0009.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0010.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0011.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0012.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0013.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0014.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0015.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0016.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0017.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0018.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0019.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0020.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0021.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0022.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0023.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0024.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0025.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0026.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0027.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0028.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0029.xsg', 1; 
        'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\IPSC\PB0002AAAA0030.xsg', 1; 
    
          };
      
      
       case '062915_PB0002_Interpatch_After Picrotoxin'
         listout = { 
             
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\After picotoxin\PB0002AAAA0151.xsg', 1; 
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\After picotoxin\PB0002AAAA0152.xsg', 1; 
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\After picotoxin\PB0002AAAA0153.xsg', 1; 
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\After picotoxin\PB0002AAAA0154.xsg', 1; 
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\After picotoxin\PB0002AAAA0155.xsg', 1; 
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\After picotoxin\PB0002AAAA0156.xsg', 1; 
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\After picotoxin\PB0002AAAA0157.xsg', 1; 
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\After picotoxin\PB0002AAAA0158.xsg', 1; 
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\After picotoxin\PB0002AAAA0159.xsg', 1; 
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\After picotoxin\PB0002AAAA0160.xsg', 1; 
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\After picotoxin\PB0002AAAA0161.xsg', 1; 
         
         };
     
      case '062915_PB0002_Interpatch_EPSC'
         listout = { 
             
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0031.xsg', 2; 
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0033.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0035.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0041.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0042.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0043.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0046.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0058.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0060.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0061.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0067.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0069.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0070.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0071.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0072.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0073.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0074.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0075.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0077.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0082.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0083.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0084.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0085.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0089.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0091.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\062915\PB0002\EPSC\PB0002AAAA0099.xsg', 2;
         };
     
     
     
    
    case '070915_PB0001_Interpatch_CPP only'
         listout = { 
             
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0011.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0012.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0013.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0014.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0015.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0016.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0017.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0018.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0019.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0020.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0021.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0022.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0023.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0024.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0025.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0026.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0027.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0028.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0029.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0030.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0031.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0032.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0033.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0034.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0035.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0036.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0037.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0038.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\PV2Pyr\PB0001AAAA0039.xsg', 1;
         
         
         
         };
         
     
     case '070915_PB0001_Interpatch_CPP and DNQX'
         listout = { 
             
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0100.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0101.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0102.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0103.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0104.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0105.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0107.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0108.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0109.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0110.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0111.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0112.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0113.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0114.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0115.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0116.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0117.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0118.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0119.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0120.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0121.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0121.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0122.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0123.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0124.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0125.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0126.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0127.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0128.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0129.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0130.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0131.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0132.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0133.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0134.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0135.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0136.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0137.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0138.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0139.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0140.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0141.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0142.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0143.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0144.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0145.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0146.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0147.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\PV2Pyr\PB0001AAAA0148.xsg', 1;
         
         };
     
     
     case '070915_PB0001_Interpatch_Picrotoxin'
         listout = { 
             
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0276.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0277.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0278.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0279.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0280.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0281.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0282.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0283.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0284.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0285.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0286.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0287.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0288.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0289.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0290.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0291.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0292.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0293.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\Picrotoxin\PB0001AAAA0294.xsg', 1;
         
         };
     
     case '070915_PB0001_Interpatch_EPSC'
         listout = { 
             
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\Pyr2PV\PB0001AAAA0040.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\Pyr2PV\PB0001AAAA0041.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\Pyr2PV\PB0001AAAA0044.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\Pyr2PV\PB0001AAAA0045.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\Pyr2PV\PB0001AAAA0049.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\Pyr2PV\PB0001AAAA0052.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\Pyr2PV\PB0001AAAA0053.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\Pyr2PV\PB0001AAAA0054.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\Pyr2PV\PB0001AAAA0055.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\only CPP\Pyr2PV\PB0001AAAA0058.xsg', 2;
         
         };
         
     case '070915_PB0001_Interpatch_EPSC_DNQX block'
         listout = { 
             
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\Pyr2PV\PB0001AAAA0230.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\Pyr2PV\PB0001AAAA0231.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\Pyr2PV\PB0001AAAA0232.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\Pyr2PV\PB0001AAAA0233.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\Pyr2PV\PB0001AAAA0234.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\070915\PB0001\CPP and DNQX\Pyr2PV\PB0001AAAA0235.xsg', 2;
          };
      
      
      
     
    case '071515_PB0001_patch_IPSC'
         listout = {  
             
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0010.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0011.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0012.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0013.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0014.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0015.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0016.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0017.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0018.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0019.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0020.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0021.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0022.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0023.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0024.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0025.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0026.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0027.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0028.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0029.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0030.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0031.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0032.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0035.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0036.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0037.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0038.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0039.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0040.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0041.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0042.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0043.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0044.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0045.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0046.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0048.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0049.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\PV2Pyr\PB0001AAAA0050.xsg', 1;
         };
     
     
       case '071515_PB0001_patch_EPSC'
         listout = {    
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\Pyr2PV\PB0001AAAA0091.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\Pyr2PV\PB0001AAAA0095.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\Pyr2PV\PB0001AAAA0096.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\Pyr2PV\PB0001AAAA0101.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\Pyr2PV\PB0001AAAA0103.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\Pyr2PV\PB0001AAAA0105.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\Pyr2PV\PB0001AAAA0106.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\Pyr2PV\PB0001AAAA0109.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\Pyr2PV\PB0001AAAA0111.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\Pyr2PV\PB0001AAAA0113.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\Pyr2PV\PB0001AAAA0115.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\Pyr2PV\PB0001AAAA0123.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\Pyr2PV\PB0001AAAA0127.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\Pyr2PV\PB0001AAAA0130.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\Pyr2PV\PB0001AAAA0138.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\Pyr2PV\PB0001AAAA0141.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0001\Pyr2PV\PB0001AAAA0142.xsg', 2;
         
         
         };
     
     case '071515_PB0003_Interpatch_IPSC'
         listout = { 
             
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\PV2Pyr\PB0003AAAA0001.xsg', 1;
         
         };
     
     
     case '071515_PB0003_Interpatch_EPSC'
         listout = { 
             
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0057.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0064.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0065.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0070.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0074.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0082.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0083.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0085.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0089.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0091.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0092.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0093.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0094.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0096.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0098.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0099.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0101.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0107.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071515\PB0003\Pyr2PV\PB0003AAAA0111.xsg', 2;
         
         
         };
     
     case '071615_PB0001_patch_IPSC'
         listout = { 
             
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0041.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0042.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0046.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0047.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0048.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0050.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0051.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0052.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0053.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0054.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0055.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0057.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0058.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0059.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0060.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0061.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0062.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0063.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0064.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0065.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0066.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0067.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0068.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0069.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\PV2Pyr\PB0001AAAA0070.xsg', 1;
         };
     
      case '071615_PB0001_patch_EPSC'
         listout = { 
             
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\Pyr2PV\PB0001AAAA0071.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\Pyr2PV\PB0001AAAA0072.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\Pyr2PV\PB0001AAAA0073.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\Pyr2PV\PB0001AAAA0074.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\Pyr2PV\PB0001AAAA0077.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\Pyr2PV\PB0001AAAA0080.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\Pyr2PV\PB0001AAAA0082.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\Pyr2PV\PB0001AAAA0084.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\Pyr2PV\PB0001AAAA0086.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\Pyr2PV\PB0001AAAA0087.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\Pyr2PV\PB0001AAAA0088.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\Pyr2PV\PB0001AAAA0091.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0001\Pyr2PV\PB0001AAAA0099.xsg', 2;
         
         };
     
     case '071615_PB0002_patch_IPSC'
         listout = { 
             
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0031.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0032.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0033.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0035.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0038.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0039.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0041.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0042.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0043.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0045.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0046.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0047.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0050.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0051.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0052.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0053.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0054.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0055.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0056.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0057.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0058.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0059.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0061.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0062.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0063.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0064.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0065.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0002\PV2Pyr\PB0002AAAA0067.xsg', 1;
         
       
         };
     
     case '071615_PB0003_Interpatch_IPSC'
         listout = { 
             
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0003\PV2Pyr\After DNQX\PB0003AAAA0116.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0003\PV2Pyr\After DNQX\PB0003AAAA0117.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0003\PV2Pyr\After DNQX\PB0003AAAA0118.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0003\PV2Pyr\After DNQX\PB0003AAAA0119.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0003\PV2Pyr\After DNQX\PB0003AAAA0131.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0003\PV2Pyr\After DNQX\PB0003AAAA0132.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0003\PV2Pyr\After DNQX\PB0003AAAA0133.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0003\PV2Pyr\After DNQX\PB0003AAAA0135.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0003\PV2Pyr\After DNQX\PB0003AAAA0136.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0003\PV2Pyr\After DNQX\PB0003AAAA0137.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0003\PV2Pyr\After DNQX\PB0003AAAA0138.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0003\PV2Pyr\After DNQX\PB0003AAAA0139.xsg', 2;
         
         
         };
     
      case '071615_PB0004_Patch_IPSC'
         listout = { 
             
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0004\PB0004AAAA0001.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0004\PB0004AAAA0002.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0004\PB0004AAAA0003.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0004\PB0004AAAA0004.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0004\PB0004AAAA0006.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0004\PB0004AAAA0008.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0004\PB0004AAAA0009.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0004\PB0004AAAA0011.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0004\PB0004AAAA0012.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0004\PB0004AAAA0014.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0004\PB0004AAAA0015.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0004\PB0004AAAA0016.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0004\PB0004AAAA0019.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071615\PB0004\PB0004AAAA0025.xsg', 1;
          };
      
      case '071715_PB0001_Patch_IPSC'
         listout = { 
         
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0001\PV2Pyr\PB0001AAAA0036.xsg', 1;
       
         };
     
     case '071715_PB0003_InterPatch_IPSC'
         listout = { 
         
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0003\PV2Pyr\PB0003AAAA0007.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0003\PV2Pyr\PB0003AAAA0048.xsg', 1;
       
         };
     
     case '071715_PB0005_Interpatch_IPSC'
         listout = { 
         
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0005\PB0005AAAA0001.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0005\PB0005AAAA0002.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0005\PB0005AAAA0004.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0005\PB0005AAAA0006.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0005\PB0005AAAA0009.xsg', 1;
       
         };
     
     case '071715_PB0006_Interpatch_IPSC'
         listout = { 
         
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0004.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0005.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0007.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0009.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0010.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0011.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0012.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0013.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0015.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0017.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0018.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0019.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0020.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0021.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0024.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0026.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0033.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0035.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0036.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0037.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0040.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0041.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0006\PV2Pyr\PB0006AAAA0043.xsg', 1;
       
         };
    
    case '071715_PB0007_Interatch_IPSC'
         listout = { 
         
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0007\PV2Pyr\PB0007AAAA0019.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0007\PV2Pyr\PB0007AAAA0020.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0007\PV2Pyr\PB0007AAAA0022.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0007\PV2Pyr\PB0007AAAA0031.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0007\PV2Pyr\PB0007AAAA0032.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0007\PV2Pyr\PB0007AAAA0046.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0007\PV2Pyr\PB0007AAAA0048.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0007\PV2Pyr\PB0007AAAA0057.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0007\PV2Pyr\PB0007AAAA0060.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0007\PV2Pyr\PB0007AAAA0064.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0007\PV2Pyr\PB0007AAAA0069.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0007\PV2Pyr\PB0007AAAA0072.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0007\PV2Pyr\PB0007AAAA0076.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0007\PV2Pyr\PB0007AAAA0077.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0007\PV2Pyr\PB0007AAAA0160.xsg', 1;
         
       
         };
     
      case '071715_PB0009_Interpatch_IPSC'
         listout = { 
         
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071715\PB0009\PB0009AAAA0001.xsg', 1;
       
         };
     
     case '071915_PB0001_Patch_IPSC'
         listout = { 
         
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0001\PV2Pyr\PB0001AAAA0012.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0001\PV2Pyr\PB0001AAAA0021.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0001\PV2Pyr\PB0001AAAA0029.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0001\PV2Pyr\PB0001AAAA0030.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0001\PV2Pyr\PB0001AAAA0057.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0001\PV2Pyr\PB0001AAAA0062.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0001\PV2Pyr\PB0001AAAA0074.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0001\PV2Pyr\PB0001AAAA0077.xsg', 1;
       
         };
     
     case '071915_PB0002_Patch_IPSC'
         listout = { 
         
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0001.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0002.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0004.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0005.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0006.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0007.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0008.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0011.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0012.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0013.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0014.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0016.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0018.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0019.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0020.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0022.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0024.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0025.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0028.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0029.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0030.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0031.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0032.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0033.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0038.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0039.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0045.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0046.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\071915\PB0002\PV2Pyr\PB0002AAAA0050.xsg', 1;
         
         };
     
    case '072215_PB0001_Interpatch_IPSC'
         listout = { 
         
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\PV2Pyr\PB0001AAAA0120.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\PV2Pyr\PB0001AAAA0121.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\PV2Pyr\PB0001AAAA0122.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\PV2Pyr\PB0001AAAA0123.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\PV2Pyr\PB0001AAAA0124.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\PV2Pyr\PB0001AAAA0125.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\PV2Pyr\PB0001AAAA0126.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\PV2Pyr\PB0001AAAA0127.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\PV2Pyr\PB0001AAAA0128.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\PV2Pyr\PB0001AAAA0129.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\PV2Pyr\PB0001AAAA0130.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\PV2Pyr\PB0001AAAA0131.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\PV2Pyr\PB0001AAAA0132.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\PV2Pyr\PB0001AAAA0133.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\PV2Pyr\PB0001AAAA0134.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\PV2Pyr\PB0001AAAA0135.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\PV2Pyr\PB0001AAAA0136.xsg', 1;
         };
     
     case '072215_PB0001_Interpatch_EPSC'
         listout = { 
         
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0020.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0021.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0023.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0025.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0028.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0029.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0032.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0034.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0036.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0039.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0061.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0062.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0063.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0065.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0068.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0069.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0071.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0075.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0078.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0084.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0086.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0098.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0103.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0110.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0119.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0211.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0217.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0001\Pyr2PV\PB0001AAAA0219.xsg', 2;
         
         };
     
     case '072215_PB0002_Interpatch_IPSC'
         listout = { 
         
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0110.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0111.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0112.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0113.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0114.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0115.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0116.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0117.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0118.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0119.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0120.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0121.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0122.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0123.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0124.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0125.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0126.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0127.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0128.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0129.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0130.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0131.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0132.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0133.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0134.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0135.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0136.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0137.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0138.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0139.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0140.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0141.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0142.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0143.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0144.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0145.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0146.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0147.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0148.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0149.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0150.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0151.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0152.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0153.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0154.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0155.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0156.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0157.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0158.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0159.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0160.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0161.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0162.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0163.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0164.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0165.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0166.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0167.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0168.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0169.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0171.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0173.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0176.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0177.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0178.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0180.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0181.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0182.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0183.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0184.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0188.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0189.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0191.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0192.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0193.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0194.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0195.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0200.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\PV2Pyr\PB0002AAAA0201.xsg', 1;
              
    
         };
     
     case '072215_PB0002_Interpatch_EPSC'
         listout = { 
         
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0010.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0011.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0012.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0020.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0030.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0031.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0034.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0035.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0047.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0048.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0049.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0051.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0057.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0080.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0107.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0210.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0218.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0220.xsg', 2;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0002\Pyr2PV\PB0002AAAA0241.xsg', 2;
         };
     
     case '072215_PB0003_Patch_IPSC'
         listout = { 
         
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0003\PV2Pyr\PB0003AAAA0001.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0003\PV2Pyr\PB0003AAAA0002.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0003\PV2Pyr\PB0003AAAA0015.xsg', 1;
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0003\PV2Pyr\PB0003AAAA0041.xsg', 1;
         };
     
      case '072215_PB0003_Patch_EPSC'
         listout = { 
         
         'C:\Ephus\R210\Ephus\Data\PawanData\Analysis\Paired recording\072215\PB0003\Pyr2PV\PB0003AAAA0051.xsg', 1;
         };
         
end




end

