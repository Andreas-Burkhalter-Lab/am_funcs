%%%% create table out of excel file, treating the first row as variablenames
% rows with empty/NaN in the first column will be deleted
%%%%% updated 18-07-02 on thermaltake

function tableout = excel2table(xlsfile)

[~,~,raw] = xlsread(xlsfile);
varnames = raw(1,:); % get varnames from first line
raw = raw(2:end,:);   % delete first line
tableout = cell2table(raw,'VariableNames',varnames);
tableout = tableout(~isnan(tableout{:,1}),:); % delete NaN rows