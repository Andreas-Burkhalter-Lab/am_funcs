function hdarray=get_high_density_array(side)
% hdarray=get_high_density_array(side)
%
% PRE:
%    side: 'left' or 'right'
% POST:
%    hdarray: a 6x5 matrix which is the high density array channels.
% NOTE:
%    the channels are numbered in Labview's way used in Berry lab.
%    For our lab, use get_hda_chan().

   hdarray=[...
       5  6  8 10 11; ...
       3  4  9 12 13; ...
       1  2  7 14 15; ...
      30 29 24 17 16; ...
      28 27 22 19 18; ...
      26 25 23 21 20; ...
      ];
   
   if(isequal(side, 'right'))
      hdarray=hdarray+30;
      hdarray=hdarray(end:-1:1,:); % flip upside down
      hdarray=hdarray(:,end:-1:1); % flip left to right
   end
   
   % -1 to get comedi channel numbering
   hdarray=hdarray-1;
   
