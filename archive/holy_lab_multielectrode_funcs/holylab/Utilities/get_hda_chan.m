function hdarray=get_hda_chan(side)
% hdarray=get_high_density_array(side)
% PRE:
%    side: 1 or 2 or 'all'
% POST:
%    hdarray: a 1x30 vector which is the high density array channels.
% NOTE:
%    the channels are numbered in comedi's way used in our lab.

   hdarray=[30 31 29 28 23 22 15 21 14 6 13 5 20 12 4 3 11 19 2 10 1 9 18 8 17 16 27 26 24 25 ...
      33 32 34 35 40 41 48 42 49 57 50 58 43 51 59 60 52 44 61 53 62 54 45 55 46 47 36 37 39 38];

   if(side==1)
      hdarray=hdarray(1:30);
   elseif(side==2)
      hdarray=hdarray(31:60)
   elseif(~isequal(side, 'all'))
      error('Unknown side');
   end
   
