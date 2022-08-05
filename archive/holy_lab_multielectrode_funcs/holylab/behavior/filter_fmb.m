function fmb_out = filter_fmb(fmb_in, origin)
%FILTER_FMB adds filtered version of fmb data to fmb struct
% Syntax:
%  fmb_out = filter_fmb(fmb_in, origin)
%     fmb_in can be a string containing the ".fmbx" file containing a fmb
%     structure containing fmb_matrix 
%     origin is a string containing the name of the PC which recorded the
%     data.  This is important for ensuring the time base is correct for
%     future plotting and analysis
%
% See also parse_fmb, nonsparse_fmb, plot_fmb

% Copyright 2008 Julian P Meeks (Timothy Holy Laboratory)

%% check input
if ischar(fmb_in)
    if ~exist(fmb_in, 'file')
        error('fmb_in is not a valid file in the current directory.');
    else
        load(fmb_in, '-mat');
        savename = fmb_in;
        save_it = true;
        if exist('fmb_matrix', 'var')
            fmb.matrix = fmb_matrix;
        end
        fmb_in = fmb;
    end
else
    save_it = false;
    fprintf('nonsparse_fmb note: data will not be saved automatically, please manually save output data in your desired file ASAP.');
end

if ~isfield(fmb_in, 'nonsparse')
    error('fmb_in does not contain the required ''nonsparse'' subfield. Please use ''nonsparse_fmb'' first.');
end

%% Loop (preliminary way to do this) to de-sparsify
if ~isfield(fmb_in, 'filtered')
    samplefreq = 100;  % sample frequency in Hz
    order = 5;
    range = 1;         % passband stop for lowpass filter in Hz
    [b, a] = cheby1(order, 0.5,...
        range/(samplefreq/2), 'low');
    
    fmb_out = fmb_in;

    for idx = 1:4
        tofilter = fmb_in.nonsparse{idx}(2,:);
        fmb_out.filtered{idx}(2,:) = filter(b,a ,tofilter);
        fmb_out.filtered{idx}(1,:) = fmb_in.nonsparse{idx}(1,:);
    end

    if save_it == true
        fmb = fmb_out;
        save(savename, 'fmb');
    end
else
    fprintf('fmb_in already contains a ''filtered'' field.  filter_fmb aborted');
    fmb_out = fmb_in;
end
