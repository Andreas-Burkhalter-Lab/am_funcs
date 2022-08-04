%boostrap to detect percent correct
%express the mean percent correct at all coherences
bootstp_num = 1000
for b = 1 : bootstp_num
    % bootstrap dataset first
    for k = 1: length(unique_coherence)
        select_boot = logical((coherence == unique_coherence(k))
        fraction correct = perc corr(select_boot); % percent correct behavioral data
        fraction_bootstrap(j) = percent correct(1); %always take the first one element
    end