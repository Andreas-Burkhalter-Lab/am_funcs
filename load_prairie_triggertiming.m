%%%% get locomotion data recorded on treadmill at 4th floor scope
% updated 2018-12-04 on thermaltake

function trigdata = load_prairie_triggertiming(trigger_data_file_name)

% import and arrange data
if strcmp(trigger_data_file_name(end-3:end),'.csv')
    trigdata = importdata(trigger_data_file_name);
    varnames = trigdata.colheaders;
    varnames = cellfun(@(x)strrep(x,'(','_'),varnames,'UniformOutput',0);
    varnames = cellfun(@(x)strrep(x,')',''),varnames,'UniformOutput',0);
    varnames = cellfun(@(x)strrep(x,' ',''),varnames,'UniformOutput',0);
    trigdata = array2table(trigdata.data,'VariableNames',varnames);
else
    fprintf(['\n Trigger timing file ' trigger_data_file_name ' not readable as .csv; opening as .mat instead.\n'])
    load(trigger_data_file_name, '-mat') %%% should have table named 'trigdata'
end
