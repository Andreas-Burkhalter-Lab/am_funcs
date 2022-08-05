function singleT = make_shifts(T,shift)
% MAKE_SHIFTS: create shifted versions of templates
% Syntax:
%   Tshift = make_shifts(T,shift)
% where
%   T is a n_channels-by-len matrix (len = length in scans of template)
%   shift is a vector containing the amount by which to shift left/right,
%     e.g., -5:5. A negative number means the template is shifted left, a
%     positive number to the right.
% and
%   Tshift is a 3d array of singles, n_channels-by-len-by-n_shifts; this is
%     the format expected by the new MultiChannelFit.
%
% You can call this for each distinct template, and concatenate. Don't
% forget to add the zero template at the end! You can do this simply by
%    Tshift(:,:,end+1) = 0;

% Copyright 2007 by Timothy E. Holy

template_length = size(T,2);
n_channels = size(T,1);
if ~isequal(shift,round(shift))
  error('For now, shift has to be an integer');
end
singleT = zeros([n_channels,template_length,length(shift)],'single');
for i = 1:length(shift)
  rng1 = make_shift_range(template_length,shift(i));
  rng2 = make_shift_range(template_length,-shift(i));
  singleT(:,rng1,i) = T(:,rng2);
end

function rng = make_shift_range(len,shift)
  rng = (1:len) + shift;
  rng = rng(rng > 0 & rng <= len);
