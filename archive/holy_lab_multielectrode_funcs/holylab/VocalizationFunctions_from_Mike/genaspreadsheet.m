function genaspreadsheet(data,cutoff,outfilename)
% Outputs Excel sheet for Gena Konopka's analyses
% Requires data cell array from pupsong_pipeline scripts

header = cat(2,data(1,:));

% Excel Spreadsheet header and layout

xlsheads = {'Id',...
    'Day',...
    'Whistle Number',...
    'Mean Duration',...
    'Min Duration',...
    'Max Duration',...
    'Mean Pause',...
    'Min Pause',...
    'Max Pause',...
    'Bout Number',...
    'Fraction of Bouts size=1',...
    'Fraction of Bouts size>1',...
    'Mean Bout Duration',...
    'Min Bout Duration',...
    'Max Bout Duration',...
    'Mean Intrabout Pause',...
    'Min Intrabout Pause',...
    'Max Intrabout Pause',...
    'Mean Interbout Pause',...
    'Min Interbout Pause',...
    'Max Interbout Pause',...
    'Fraction of SS Whistles',...
    'Mean SS Whistle Mean frequency',...
    'Min SS Whistle Mean frequency',...
    'Max SS Whistle Mean frequency',...
    'Mean SS Whistle Median frequency',...
    'Min SS Whistle Median frequency',...
    'Max SS Whistle Median frequency',...
    'Mean SS Whistle Frequency Range',...
    'Min SS Whistle Frequency Range',...
    'Max SS Whistle Frequency Range',...
    'Mean SS Whistle Slope',...
    'Min SS Whistle Slope',...
    'Max SS Whistle Slope',...
    'Fraction of Jump Whistles',...
    'Mean Jump Whistle Mean frequency',...
    'Min Jump Whistle Mean frequency',...
    'Max Jump Whistle Mean frequency',...
    'Mean Jump Whistle Median frequency',...
    'Min Jump Whistle Median frequency',...
    'Max Jump Whistle Median frequency',...
    'Mean Jump Whistle Frequency Range',...
    'Min Jump Whistle Frequency Range',...
    'Max Jump Whistle Frequency Range',...
    'Jump Whistle - Mean Number of Jumps',...
    'Jump Whistle - Min Number of Jumps',...
    'Jump Whistle - Max Number of Jumps'};

ncols = size(xlsheads,2);
nrows = size(data,1) - 1;
headidx = zeros(1,ncols);
outputTable = zeros(nrows,ncols);
    
headidx(1) = find(strcmp(header,'id'));
headidx(2) = find(strcmp(header,'day'));
headidx(3) = find(strcmp(header,'whisn'));
headidx(4:6) = find(strcmp(header,'dt'));
headidx(7:9) = find(strcmp(header,'pause'));
headidx(10) = find(strcmp(header,'boutn'));
headidx(11:12) = find(strcmp(header,'boutsize'));
headidx(13:15) = find(strcmp(header,'boutdt'));
headidx(16:21) = find(strcmp(header,'pause'));
headidx(22) = find(strcmp(header,'f(ss)'));
headidx(23:25) = find(strcmp(header,'mf(ss)'));
headidx(26:28) = find(strcmp(header,'medf(ss)'));
headidx(29:31) = find(strcmp(header,'fr(ss)'));
headidx(32:34) = find(strcmp(header,'slope(ss)'));
headidx(35) = find(strcmp(header,'f(jumps)'));
headidx(36:38) = find(strcmp(header,'mf(jumps)'));
headidx(39:41) = find(strcmp(header,'medf(jumps)'));
headidx(42:44) = find(strcmp(header,'fr(jumps)'));
headidx(45:47) = find(strcmp(header,'jumpn'));

    for i = 1:nrows
        for j = 1:ncols
        x = data{i+1,headidx(j)};
            if any(j == [1,2,3,10,22,35])
                outputTable(i,j) = x;
            elseif any(j == [4,7,13,23,26,29,32,36,39,42,45])
                outputTable(i,j) = mean(x);
            elseif any(j == [5,8,14,24,27,30,33,37,40,43,46])
                if ~isempty(x)
                outputTable(i,j) = min(x);
                else
                outputTable(i,j) = NaN;
                end
            elseif any(j == [6,9,15,25,28,31,34,38,41,44,47])
                if ~isempty(x)
                outputTable(i,j) = max(x);
                else
                outputTable(i,j) = NaN;
                end
            elseif j == 11
                outputTable(i,j) = length(find(x==1));
                outputTable(i,j) = outputTable(i,j)/length(x);
            elseif j == 12
                outputTable(i,j) = length(find(x>1));
                outputTable(i,j) = outputTable(i,j)/length(x);
            elseif any(j == [16,17,18,19,20,21])
                if data{i+1,2} == 4
                    if j == 16
                        x = x(x<=cutoff.day4);
                        outputTable(i,j) = mean(x);
                    elseif j == 17
                        x = x(x<=cutoff.day4);
                        if ~isempty(x)
                        outputTable(i,j) = min(x);
                        else
                        outputTable(i,j) = NaN;
                        end
                    elseif j == 18
                        x = x(x<=cutoff.day4);
                        if ~isempty(x)
                        outputTable(i,j) = max(x);
                        else
                        outputTable(i,j) = NaN;
                        end
                    elseif j == 19
                        x = x(x>cutoff.day4);
                        outputTable(i,j) = mean(x);
                    elseif j == 20
                        x = x(x>cutoff.day4);
                        if ~isempty(x)
                        outputTable(i,j) = min(x);
                        else
                        outputTable(i,j) = NaN;
                        end
                    elseif j == 21
                        x = x(x>cutoff.day4);
                        if ~isempty(x)
                        outputTable(i,j) = max(x);
                        else
                        outputTable(i,j) = NaN;
                        end
                    end                        
                elseif data{i+1,2} == 7
                    if j == 16
                        x = x(x<=cutoff.day7);
                        outputTable(i,j) = mean(x);
                    elseif j == 17
                        x = x(x<=cutoff.day7);
                        if ~isempty(x)
                        outputTable(i,j) = min(x);
                        else
                        outputTable(i,j) = NaN;
                        end
                    elseif j == 18
                        x = x(x<=cutoff.day7);
                        if ~isempty(x)
                        outputTable(i,j) = max(x);
                        else
                        outputTable(i,j) = NaN;
                        end
                    elseif j == 19
                        x = x(x>cutoff.day7);
                        outputTable(i,j) = mean(x);
                    elseif j == 20
                        x = x(x>cutoff.day7);
                        if ~isempty(x)
                        outputTable(i,j) = min(x);
                        else
                        outputTable(i,j) = NaN;
                        end
                    elseif j == 21
                        x = x(x>cutoff.day7);
                        if ~isempty(x)
                        outputTable(i,j) = max(x);
                        else
                        outputTable(i,j) = NaN;
                        end
                    end
                elseif data{i+1,2} == 10
                    if j == 16
                        x = x(x<=cutoff.day10);
                        outputTable(i,j) = mean(x);
                    elseif j == 17
                        x = x(x<=cutoff.day10);
                        if ~isempty(x)
                        outputTable(i,j) = min(x);
                        else
                        outputTable(i,j) = NaN;
                        end
                    elseif j == 18
                        x = x(x<=cutoff.day10);
                        if ~isempty(x)
                        outputTable(i,j) = max(x);
                        else
                        outputTable(i,j) = NaN;
                        end
                    elseif j == 19
                        x = x(x>cutoff.day10);
                        outputTable(i,j) = mean(x);
                    elseif j == 20
                        x = x(x>cutoff.day10);
                        if ~isempty(x)
                        outputTable(i,j) = min(x);
                        else
                        outputTable(i,j) = NaN;
                        end
                    elseif j == 21
                        x = x(x>cutoff.day10);
                        if ~isempty(x)
                        outputTable(i,j) = max(x);
                        else
                        outputTable(i,j) = NaN;
                        end
                    end
                elseif data{i+1,2} == 14
                    if j == 16
                        x = x(x<=cutoff.day14);
                        outputTable(i,j) = mean(x);
                    elseif j == 17
                        x = x(x<=cutoff.day14);
                        if ~isempty(x)
                        outputTable(i,j) = min(x);
                        else
                        outputTable(i,j) = NaN;
                        end
                    elseif j == 18
                        x = x(x<=cutoff.day14);
                        if ~isempty(x)
                        outputTable(i,j) = max(x);
                        else
                        outputTable(i,j) = NaN;
                        end
                    elseif j == 19
                        x = x(x>cutoff.day14);
                        outputTable(i,j) = mean(x);
                    elseif j == 20
                        x = x(x>cutoff.day14);
                        if ~isempty(x)
                        outputTable(i,j) = min(x);
                        else
                        outputTable(i,j) = NaN;
                        end
                    elseif j == 21
                        x = x(x>cutoff.day14);
                        if ~isempty(x)
                        outputTable(i,j) = max(x);
                        else
                        outputTable(i,j) = NaN;
                        end
                    end
                end
            end
        end
    end

    
    % Write to excel spreadsheet
    
    xlswrite(outfilename,outputTable,'output','A2');
    xlswrite(outfilename,xlsheads,'output','A1');
end
