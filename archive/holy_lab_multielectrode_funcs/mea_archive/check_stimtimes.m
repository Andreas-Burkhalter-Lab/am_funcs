function [ varargout ] = check_stimtimes( snipstruct_all)
%CHECK_STIMTIMES Graphically check for correct stim time assignments. 
%   Requires a 'snipstruct_all' structure as output by see_pre_post,
%   including 'stimwave' fields as created by spike_vs_stim

rows = 5; 
cols = 5;

for i = 1:length(snipstruct_all)
    subplot(rows,cols,i)
    plot(snipstruct_all(i).stimwave,'b')
    hold('on');
    scatter(snipstruct_all(i).stimtime,snipstruct_all(i).stimthresh,'r')
end

end

