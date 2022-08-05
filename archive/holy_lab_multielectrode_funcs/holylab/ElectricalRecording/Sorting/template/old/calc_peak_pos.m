function peakPos=calc_peak_pos(spike, sniprange)
% this one-liner is useful when we do have off-by-one problem
% TODO: replace the relevent code in verify_fit()
% TODO: add comment

   peakPos=spike.shiftedTime;
   