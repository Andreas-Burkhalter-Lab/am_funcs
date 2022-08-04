function [ varargout ] = see_prepost( snipstruct_all )
%SEE_PREPOST Plot response timing info. 
%   Requires snipstruct_all variable as created by collect_snipstruct.

rows = 5;
cols = 5;

for i = 1:length(snipstruct_all)    %% i = each file or element of snipstruct_all
    subplot(rows,cols,i)
    hold('on');
    for j = 1:length(snipstruct_all(i).channels)    %% j = each channel
        scatter(snipstruct_all(i).tpre{j},snipstruct_all(i).channels(j)*ones(1,length(snipstruct_all(i).tpre{j})),'b')
        scatter(snipstruct_all(i).tpost{j},snipstruct_all(i).channels(j)*ones(1,length(snipstruct_all(i).tpost{j})),'r')
    end
    title(snipstruct_all(i).fh.filename(1:end-5),'Interpreter','none','FontSize',9)
    ylabel('Channel Number')
end

