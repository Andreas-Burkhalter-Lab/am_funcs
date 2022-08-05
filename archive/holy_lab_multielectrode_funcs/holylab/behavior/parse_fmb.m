function fmb = parse_fmb(filename)
% PARSE_FMB turns a .fmb file into a .fmbx file (.mat file)
% ".fmb" files are text outputs of the jstest --event function, and
%        contain lines with times, identities, and values from USB
%        joystick axes attached to "fake mouse butts" (FMBs).
%
% This program takes a file and creates a ".mat" file of the same name as
% the input filename but the suffix is ".fmbx".
% The function returns fmb.matrix, a matrix containing a cell array of
% vectors, each cell contains the vector data for a single joystick axis.

% Copyright 2008 Julian P Meeks (Timothy Holy Laboratory)
% 
% Revision History:
% 2008_05_03: Wrote it (JPM) with support for only a single joystick and 4
%             axes
% 2008_07_14: Updated script for multiple acquisition PCs

if ~exist(filename, 'file')
    error('%s does not exist.',filename);
end

in_fid = fopen(filename, 'r');  % open in read-only mode
line_count = 1;
tline = 0;


for idx = 1:30
    tline = fgets(in_fid);
end

% set up cell array to handle data:
temp_fmb_matrix = cell(1,4);

while tline ~= -1
    tline  = fgets(in_fid);
    % find axis number (and use to assign cell)
    ax_mark = strfind(tline,'number')+7;
    thisax_s = [];
    if ~isempty(ax_mark) && ax_mark <= size(tline,2)
        while ~isempty(str2num(tline(ax_mark)))
            if isempty(thisax_s)
                thisax_s = tline(ax_mark);
            else
                thisax_s(end+1) = tline(ax_mark);
            end
            if ax_mark+1 <= size(tline,2)
                ax_mark = ax_mark + 1;
            else
                thisax_s = [];
                break;
            end
        end
        if ~isempty(thisax_s)
            thisax = str2num(thisax_s);
        else
            break;
        end
        currentcell = [];
        switch thisax
            case 0
                currentcell = 1;
            case 1
                currentcell = 2;
            case 2
                currentcell = 3;
            case 3
                currentcell = 4;
        end
    
        if ~isempty(currentcell)
            % grab the timepoint number
            t_mark = strfind(tline, 'time') + 5;
            timept_s = [];
            while ~isempty(str2num(tline(t_mark)))||tline(t_mark)=='-'
                if isempty(timept_s)
                    timept_s = tline(t_mark);
                else
                    timept_s(end+1) = tline(t_mark);
                end
                if t_mark+1 <= size(tline,2)
                    t_mark = t_mark + 1;
                else
                    timept_s = [];
                    break;
                end
            end
       
            % now grab the value of the axis
            v_mark = strfind(tline, 'value') + 6;
            value_s = [];
            if v_mark <= size(tline, 2)
                while ~isempty(str2num(tline(v_mark)))||tline(v_mark)=='-'
                    if isempty(value_s)
                        value_s = tline(v_mark);
                    else
                        value_s(end+1) = tline(v_mark);
                    end
                    if v_mark+1 <= size(tline,2)
                        v_mark = v_mark + 1;
                    else
                        value_s = [];
                        break;
                    end
                end
            end
            if ~isempty(timept_s)&&~isempty(value_s)
                temp_fmb_matrix{1,currentcell}(1,end+1) = str2num(timept_s);
                temp_fmb_matrix{1,currentcell}(2,end) = str2num(value_s);
            else
                break;
            end
            % if we've written data,
            line_count = line_count+1;
        end
        
    end
end

fmb.matrix = temp_fmb_matrix;

% write file
basename = filename(1:strfind(filename,'.'));
writefile = [basename 'fmbx'];
save(writefile,'fmb','-mat')

end
