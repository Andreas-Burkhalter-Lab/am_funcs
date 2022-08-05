function result=get_hda_holylab(side)
% this func returns comedi channel numbers in a high density array used in our lab
% SYNTAX:
%   result=get_hda_holylab(side)
% PRE: 
%   side: 'left' or 'right'  %NB: add more comments
% POST:
%   result: a 6x5 matrix which is the high density array channels.
% NOTE:
%   the matrix is organized according to the channels' physical locations
   
  allowed_sides = {'left','right','miprobe4x4'};
  if isempty(intersect(side,allowed_sides))
    error(['Array type ' side ' not recognized']);
  end

  if strcmp(side,'miprobe4x4')
    result = [...
      47, 38, 44, 12; ...
      39,  7, 13,  4; ...
      15,  6,  5, 37; ...
      14, 46, 45, 36; ...
      ];
    return;
  end

  idx=get_high_density_array(side);
   ch=get_hda_chan('all');
   result=ch(idx+1);
   