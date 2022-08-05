function hax = plot_fmb(fmb, origin, varargin)
% PLOT_FMB takes a fmb struct (containing cell arrays) and plots a figure, returning hax
%   fmb is fmb structure with at minimum the field 'matrix'
%          --> if 'nonsparse' field is present, this will be used as it is prettier
%   origin is a string containing the name of the machine used for acquisition
%          *** this is very important for properly displaying times ***
%          valid strings (as of 2008_06_03) are:
%            -- 'opium'
%            -- 'glow'
%   varargin can contain a string for "mode" for plotting
%   - 'combine_axes' collapses by adding axes1&2 and 3&4
%   - 'deriv' plots the derivative of position (default)
%   - 'bin_deriv' plots a sliding window integration of the derivative
%

% Copyright 2008 Julian P Meeks (Timothy Holy Laboratory)
%

%%  check variables
if ischar(fmb)
    if ~exist(fmb, 'file')
        error('fmb_in is not a valid file in the current directory.');
    else
        load(fmb, '-mat');
        if exist('fmb_matrix', 'var')
            fmb.matrix = fmb_matrix;
        end
    end
end

if nargin > 1
    for idx = 1:size(varargin,2)
        if ~ischar(varargin{idx})
            error('varargin must be a string.  ''raw'', ''deriv'', ''int_deriv'', and ''combine_axes'' are only options at this time.');
        else
            this_arg = varargin{idx};
            switch this_arg
                case {'raw', 'int_deriv', 'deriv'}
                    mode = varargin{idx};
                case 'combine_axes'
                    combine = true;
            end
        end
    end
end
                        
if ~exist('mode', 'var')
    mode = 'raw';
end
if ~exist('combine', 'var')
    combine = false;
end
    
%% Get to plotting

toffset = fmb.matrix{1,1}(1,1);

% find the last event logged over the 4 controllers and assign that time to
% 'tend'
for idx = 1:4
    tend(idx) = fmb.matrix{1,idx}(1,end);
end
tend = max(tend);

% convert the data, which are time, value pairs, to continuous values with
% a consistent time base.  Values between events remain the value prior to
% the registered event
fmb_copy = cell(1,4);

if strcmp(origin, 'glow')
    tbase = 0:10:(tend-toffset); % tbase is a vector of controller-centric time base values starting at zero(8000/sec)
    tvals = toffset:10:tend;     % tvals are the values that match the time element of fmb.matrix
elseif strcmp(origin, 'opium')
    tbase = 0:8:(tend-toffset);
    tvals = toffset:8:tend;
else
    error('''plot_fmb'' does not recognize %s as a valid machine origin.', origin);
end

if ~isfield(fmb, 'nonsparse')
    fmb = nonsparse_fmb(fmb);
end

figure
hold on;
xlims = [0 0];
ylims = [-5 5];

if strcmp(origin, 'glow')
    toseconds = 1000;
elseif strcmp(origin, 'opium')
    toseconds = 800;
end


if combine == true
    nax = 2;
else
    nax = 4;
end

for idx = 1:nax
    voffset = fmb.matrix{1,idx}(2,1);
    subplot(nax,1,idx);
    hax(idx) = gca;
%     if isfield(fmb, 'filtered')
%         toplot = fmb.filtered;
%     else
if isfield(fmb, 'nonsparse')
  %      toplot = fmb.nonsparse;
 %   else
        toplot = fmb.matrix;
    end
    
    if combine == true
        temp = cell(1);
        temp{1}(1,:) = toplot{1,1}(1,:);
        temp{1}(2,:) = toplot{1,1}(2,:) + toplot{1,2}(2,:);
        temp{2}(1,:) = temp{1}(1,:);
        temp{2}(2,:) = toplot{1,3}(2,:) + toplot{1,4}(2,:);
        toplot = temp;
    end
    
    switch mode
        case 'raw'
            plot(toplot{1,idx}(1,1:end)/toseconds/3600, (toplot{1,idx}(2,:)-voffset)/337)
        case 'deriv'
            plot(toplot{1,idx}(1,2:end)/toseconds/3600,diff(((toplot{1,idx}(2,:))-voffset)/337));
        case 'int_deriv'
            count = 1;
            sl_window = zeros(1, ceil((size(toplot{1,idx}(2,:),2))/1000));
            for i_idx = 1001:1000:size(toplot{1,idx},2)-1001
                sl_window(count) = sum(abs(diff(toplot{1,idx}(2,i_idx-1000:i_idx+999))));
                count = count+1;
            end
            plot((toplot{1,idx}(1,1:1000:end))/toseconds/3600,sl_window/337);
    end     
    cur_xlim = get(gca,'xlim');
    if cur_xlim(2) > xlims(2)
        xlims = cur_xlim;
    end
    cur_ylim = get(gca,'ylim');
    if cur_ylim(2) > ylims(2) || cur_ylim(1)<ylims(1)
        ylims = cur_ylim;
    end
end

for idx = 1:nax
    set(hax(idx),'xlim',xlims);
    set(hax(idx),'ylim',ylims);
end

end