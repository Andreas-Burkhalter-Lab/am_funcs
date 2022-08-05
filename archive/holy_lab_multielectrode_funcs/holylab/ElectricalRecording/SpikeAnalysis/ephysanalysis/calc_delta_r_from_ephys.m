function delta_r = calc_delta_r_from_ephys(ephysin,valvename,delta_r_range, whichtimes, varargin)
% function calc_delta_r_from_ephys calculates firing rate changes from ephys data
%
% delta_r = calc_delta_r_from_ephys(ephysin,valvename,delta_r_range,varargin)
%    ephysin should be given after creating a structure array of ephys
%         structures by interval (using ephyssubrange, for example)
%    valvename is a character array with the exact valve label you want to
%         use to (TO DO: include an 'all' option)
%    whichtimes is a string, either 'sniptimes' or 'celltimes'
%    delta_r_range is a two-vector (in seconds) of the time you wish to
%         analyze following stimulus onset
%    delta_r will be returned as a cell containing a double array with the same number of
%         elements as intervals with your chosen valvename
%    ** if 'celltimes' is chosen, the returned value will be a cell array
%         of doubles with number of cells equaling the number of cells,
%         each with the same number of intervals as in original
%
% See also EPHYSFROMAI, EPHYSSUBRANGE,

% Copyright Julian P. Meeks 2007 (Timothy Holy Laboratory)
% History:
%    2007_11_29: (JPM) wrote it

%% Start by checking input variables for errors
if length(ephysin)<=1
    warning('This routine is currently only intended for ephys structure arrays (by interval)');
end

if isempty(valvename)
   delta_r = [];
   return;
elseif ~isstr(valvename)
   error('valvename must be in string format');
end

if length(delta_r_range)~=2
    error('delta_r_range must be a two-vector');
elseif ~isbetween(delta_r_range(1),...
        [ephysin(1).toffset...
         (ephysin(1).scanrange(2)-ephysin(1).scanrange(1)-ephysin(1).toffset)/ephysin(1).scanrate])...
   || ~isbetween(delta_r_range(2),...
        [ephysin(1).toffset...
         (ephysin(1).scanrange(2)-ephysin(1).scanrange(1)-ephysin(1).toffset)/ephysin(1).scanrate])
     error('delta_r_range must be within the interval range');
end

if nargin>4
    a = 0 % not implemented yet
end
    
%% Next, calculate the delta r
% trying to eliminate erroneous valvelabels (careful)
for idx = 1:length(ephysin(1).valvelabels)
    if isempty(ephysin(1).valvelabels{idx})
        valvelabels = ephysin(1).valvelabels;
        valvelabels(idx) = {'null'};
        [ephysin.valvelabels]=deal(valvelabels);
    end
end

valvelabels = ephysin(1).valvelabels;

for idx = 1:length(ephysin)
    if strmatch(valvename,ephysin(idx).tag, 'exact');
        if ~exist('thisdata','var')          
            thisdata = ephysin(idx);
        else
            thisdata(end+1)=ephysin(idx);
        end
    end
end

if ~exist('thisdata','var')
    delta_r = [];
    return
end

for idx = 1:size(thisdata,2);
    delta_r_scanrange{idx} = (delta_r_range-thisdata(idx).toffset)*thisdata(idx).scanrate+thisdata(idx).scanrange(1);
    base_scanrange{idx} = [thisdata(idx).scanrange(1) 
                               thisdata(idx).scanrange(1)+...
                                     (delta_r_range(2)-delta_r_range(1))*thisdata(idx).scanrate];
    if whichtimes == 'sniptimes'
        numspikes(idx) = length(between(thisdata(idx).sniptimes{:},...
                             delta_r_scanrange{idx}(1),...
                             delta_r_scanrange{idx}(2)));
        basespikes(idx) = length(between(thisdata(idx).sniptimes{:},...
                             base_scanrange{idx}(1),...
                             base_scanrange{idx}(2)));
        stim_firerate = numspikes(idx)/(delta_r_range(2)-delta_r_range(1));
        base_firerate = basespikes(idx)/(delta_r_range(2)-delta_r_range(1));
        delta_r(idx) = stim_firerate-base_firerate;

    elseif whichtimes == 'celltimes'
        for cell_idx = 1:size(ephysin(1).cellnums,2)
            numspikes(idx) = length(between(thisdata(idx).celltimes{cell_idx},...
                delta_r_scanrange{idx}(1),...
                delta_r_scanrange{idx}(2)));
            basespikes(idx) = length(between(thisdata(idx).celltimes{cell_idx},...
                base_scanrange{idx}(1),...
                base_scanrange{idx}(2)));
            cells(cell_idx).numspikes = numspikes;
            cells(cell_idx).basespikes = basespikes;
            cells(cell_idx).stim_firerate = cells(cell_idx).numspikes(idx)/(delta_r_range(2)-delta_r_range(1));
            cells(cell_idx).base_firerate = cells(cell_idx).basespikes(idx)/(delta_r_range(2)-delta_r_range(1));
            delta_r{cell_idx}(idx) = cells(cell_idx).stim_firerate-cells(cell_idx).base_firerate;
        end
    end        
end

if ~iscell(delta_r)
    delta_r = {delta_r};
end

end