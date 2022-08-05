function fmb_out = nonsparse_fmb(fmb_in, origin)
%NONSPARSE FMB adds nonsparse version of fmb data to fmb struct
% Syntax:
%  fmb_out = nonsparse_fmb(fmb_in, origin)
%     fmb_in can be a string containing the ".fmbx" file containing a fmb
%     structure containing fmb_matrix 
%     origin is a string containing the name of the PC which recorded the
%     data.  This is important for ensuring the time base is correct for
%     future plotting and analysis
%
% See also parse_fmb, plot_fmb

% Currently uses lengthy loops to assign data values, can probably be
% optimized
%
% Copyright 2008 Julian P Meeks (Timothy Holy Laboratory)

%% check input
if iscell(fmb_in)
    max_in = size(fmb_in,2);
else
    fmb_in = {fmb_in};
    max_in = 1;
end
for idx = 1:max_in  
    if ischar(fmb_in{idx})
        if ~exist(fmb_in{idx}, 'file')
            error('%s is not a valid file in the current directory.', fmb_in{idx});
        else
            a(idx) = load(fmb_in{idx}, '-mat');
            savename{idx} = fmb_in{idx};
            savename{idx}(end) = 'n';
            save_it = true;
            if exist('fmb_matrix', 'var')
                a(idx).fmb.matrix = a(idx).fmb_matrix;
            end
            fmb_in{idx} = a(idx).fmb;
        end
    else
        save_it = false;
        fprintf('nonsparse_fmb note: data will not be saved automatically, please manually save output data in your desired file ASAP.');
    end
end

for idx = 1:max_in
    if ~isfield(fmb_in{idx}, 'matrix')
        error('fmb_in does not contain the required ''matrix'' subfield.');
    end
end

%% Loop (preliminary way to do this) to de-sparsify
% take the raw inputs and place them in a single cell array
fmb = struct;
count = 1;
for c_idx = 1:max_in
    for ax_idx = 1:size(fmb_in{c_idx}.matrix, 2)
        fmb.times{count} = fmb_in{c_idx}.matrix{ax_idx}(1,:);
        tstart(count) = fmb.times{count}(1);
        tend(count) = fmb.times{count}(end);
        fmb.events{count} = fmb_in{c_idx}.matrix{ax_idx}(2,:);
        estart(count) = fmb.events{count}(1);
        emax(count) = max(fmb.events{count});
        emin(count) = min(fmb.events{count});
        count = count+1;
    end
end
if size(tstart,2) == size(estart,2)
    nchans = size(tstart,2);
else
    error('the number of channels is not consistent.  check input parsed files for errors.');
end

for idx = 1:nchans
    if tstart(idx) ~= min(tstart);
        fmb.times{idx}(2:end+1) = fmb.times{idx}(1:end);
        fmb.times{idx}(1) = min(tstart);
        fmb.events{idx}(2:end+1) = fmb.events{idx}(1:end);
        fmb.events{idx}(1) = estart(idx);
    end
end

if strcmp(origin, 'glow')
    tbase = 0:10:(max(tend)-min(tstart)); % tbase is a vector of controller-centric time base values starting at zero(8000/sec)
    tvals = min(tstart):10:max(tend);     % tvals are the values that match the time element of fmb.matrix
elseif strcmp(origin, 'opium')
    tbase = 0:8:(max(tend)-min(tstart));
    tvals = min(tstart):8:max(tend);
else
    error('''nonsparse_fmb'' does not recognize %s as a valid machine origin.', origin);
end

data = zeros(nchans, size(tbase,2));
for idx = 1:nchans
    tic;
    for t_idx = 1:size(tbase,2)
        [c ia ib] = intersect(tvals(t_idx), fmb.times{idx});
        if ~isempty(c)
            data(idx,t_idx) = fmb.events{idx}(ib);
            current = data(idx,t_idx);
        else
            data(idx, t_idx) = current;
        end
    end
    fprintf('Axis %d: ', idx);   % can be commented
    toc;                         % can be commented
end    

% collapse the axes
if mod(nchans, 2) == 0
    combined_data = zeros(nchans/2, size(tbase,2));
    for idx = 1:nchans/2
        combined_data(idx,:) = data(1+(idx-1)*2,:)+data(2+(idx-1)*2,:);
    end
else
    fprintf('Odd number of axes, no collapsing will be saved');
end
            
fmb.nonsparse_times = tbase;
fmb.nonsparse_events = data;
fmb.combined_nonsparse_events = combined_data;

if save_it == true
    save(savename{1}, 'fmb');
    fprintf('%s saved.\n', savename{1});
end

fmb_out = fmb;