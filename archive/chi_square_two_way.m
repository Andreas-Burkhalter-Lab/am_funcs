%%%% chi-square independence test for frequency data along two dimensions
% creates table with entry for every observation so matlab can handle
%       it with crosstab.m
%%% [tableout, chi2, pval, labelsout] = chi_square_independence(freq_table, rowlabels, collabels)
%%%%% updated 18/6/14 on thermaltake

function [tableout, chi2, pval, labelsout] = chi_square_two_way(freq_table, rowlabels, collabels)

[nrows, ncols] = size(freq_table);
rowlabels = vardefault('rowlabels',cellstr(num2str([1:nrows]')));
collabels = vardefault('collabels',cellstr(num2str([-1:-ncols]')));
rowobs = cell(0); % observations corresponding to rows in freq_table
colobs = cell(0); % observations corresponding to columns in freq_table
for irow = 1:nrows
    rlab = rowlabels{irow};
    for icol = 1:ncols
        clab = collabels{icol};
        rowobs = [rowobs; repmat({rlab}, freq_table(irow,icol), 1)];
        colobs = [colobs; repmat({clab}, freq_table(irow,icol), 1)];
    end
end
rowobs = categorical(rowobs);
colobs = categorical(colobs);

[tableout, chi2, pval, labelsout] = crosstab(rowobs, colobs); % do chi square test