function realstim=timestim(filename, options)
% TIMESTIM: analyze stimulus fire times
% syntax:
%    realstim=timestim(filename, options)
% pre:
%    filename: stimulus file name
%    options: a structure with fields ---
%       wndL=0.1:  when load raw data to detect valve transition, how much data
%                 (in seconds) is read before intended transition time
%       wndR=0.1:  ... after intended transition time
%       feedback=0 or 63: the feedback channel # (comedi channel #, not
%                         virtual scsi#); default is 0 for
%                         chanel.wustl.edu, 63 for others
%       calc_vs_intended_tolerance=0-1: determines the fraction of valve
%                                   transitions for which the calcualted
%                                   valve number (based on the feedback
%                                   voltage trace) can be different from
%                                   the intendend number (based on the
%                                   header information) without timestim
%                                   choking.  NOTE: in case where there is
%                                   a discrepancy, even if it's within the
%                                   tolerance, a plot of the discrepancy
%                                   versus trial number will be generated.
%       no_intended_times = 0 or 1: If set to 1, will ignore any
%                   information from the header about what the intended
%                   stimulus sequence is, and instead attempt to parase it
%                   out entirely from the efferent copy.
%       voltage_to_vlvnum_assignments: matrix OR 'return_voltages'
%                   This will only be looked at at all if no_intended_times
%                   is set to 1.
%                   * If set to 'return_voltages' timestim won't even make
%                   an attempt to convert the voltages it comes up with
%                   into valvenumbers, and will just report the voltages
%                   directly
%                   * If you want to assign a mapping of voltages to
%                   valvenumbers, you can do it two different ways:
%                   (1) you can present the information in a two
%                   column form, with the first column indicating valve
%                   number and the second voltage;
%                   any given voltage will be assigned to whichever
%                   valve-number's voltage it's closest to, giving an
%                   error only if it's in the middle 50% of the voltage
%                   range between them.
%                   (2) you can present the information in a three column form, where
%                   the second column is the minimum voltage that
%                   corresponds to that valve number and the third is the
%                   maximum
%       testing_mode = 0 or 1: If set to 1, will output ___
%       total_override: if a "realstim" file (eg, one you pieced together
%                   by hand, painfully, with sweat and tears...) is entered
%                   here, it will be used as the output and the entire stim
%                   channel itself will be completely ignored. 
%  (see also the two more whose comments are located where the deafults are
%  set)
%
% post:
%    realstim: an nTransitions-by-2 matrix of the form [vlv time], where the column
%      vector vlv gives the valve numbers, and time gives the actual scan
%      numbers of the transitions.

% history:
%    4/26/2004: 1. change the way to "boldly guess the current valve number"
%               2. change the way to estimate the left edge of valve opening.
%    9/29/2003: handle the situation when file is shorter than stimulus
%               sequence
%    2004-07-22:  1. Adjust the expected time of a transition
%    (TEH)           incrementally as you go along, so that drift between
%                    the clocks doesn't accumulate
%                 2. Prevent problems with "fake transitions" (where the
%                    valve doesn't actually change state).
%                 3. Eliminate the need for voltage->valvenum calibration
%                 4. LFS support for large files
%                 5. Make sure consider only one transition at a time
%                    (important for transitions narrower than [-wndL wndR])
%                 6. Because of higher accuracy we can decrease the window
%                    size, making the algorithm so fast that I've removed the
%                    progress bar.
%    2006-08-10: Added in an option where you can increase the tolerance
%    (RCH)       discrepancies between the calculated and intended valve
%                numbers.  The default tolerance will remain zero; this
%                option is being designed to fit the situation where the
%                digital to analogue chip in the dumbbox seems to be
%                initially very noisy, but "settles down" shortly into the
%                experiment, such that most of the valve transitions are
%                detected correctly, providing a sufficient check that the
%                equipment is working as expected.
%    2006-08-14: Added in basically a whole new "mode" of running timestim.
%    (RCH)       Designed to be the "robust, semi-manuel" mode - so that
%                if your header file is corrupted or gone, or you don't
%                trust your intended valve times/numbers for some reason,
%                you can brute force your way through figuring out what
%                stimuli your prep actually saw. Rest of timestim should
%                (!!) work just the same, as long as
%                options.no_intended_times is left to its deafult 0...
%   2006-10-27:  Added in the total_override option; only point of putting
%    (RCH)       it here is to let you run autosnip2 with nothing but
%                timestim having to know anything about where the stimulus
%                information is coming from
%   2008-01-16:  Added in option (on by default) to allow a check anytime a
%    (RCH)       deadloop is reached to see if it's a one- or two-point
%                blip, and, if so, just remove those points...

if(nargin ~=1 & nargin~=2)
    error('timestim() requires 1 or 2 input arguments');
end

if(nargin==1)
    options=[];
end

% fill unsupplied options:
if(~isfield(options, 'wndL'))  options.wndL=0.1; end % 0.1 sec from intended position
if(~isfield(options, 'wndR'))  options.wndR=0.1; end % 0.1s from intended position
options = default(options,'no_intended_times',0); % for use in "manuel override" situation
options = default(options,'voltage_to_vlvnum_assignments','guess_mapping'); % will only be used if no_intended_times is 1;
options = default(options,'min_time_between_valve_changes',1); % time in seconds; will only be used if no_intended_times is 1
options = default(options,'calc_vs_intended_tolerance',0);
options = default(options,'nToAverage',100);
    % Set lower than the number of scans for which you
    % would actually WANT to see a fluctuation in
    % voltage
options = default(options,'testing_mode',0);
options = default(options,'total_override',0);
options = default(options,'fix_one_point_causes_of_deadloops_on',1);

if length(options.total_override) > 1
    realstim = options.total_override;
    return
end

% Open file with potential LFS support
[h,fid] = readheader(filename);

% If in testing mode, initialize a structure that will hold whatever we
% want to dump to a testing-output file at the end:
if options.testing_mode
    testing_out = struct;
    testing_out.resampled_freq_goal = 10;  % In Hz, how much resolution you want in the resampled, outputted testing data
    testing_out.resampling_factor = round(h.scanrate/testing_out.resampled_freq_goal);
    testing_out.resampling_freq = h.scanrate/testing_out.resampling_factor;
end

% check if the stimuli were deliveried by a robot
if(is_robot_stim(h))
    if(should_use_lfs(filename))
        closelfs(fid);
    else
        fclose(fid);
    end
    realstim=time_robot_stim(filename, options);
    return;
end

filemode = '';
is_force_use_no_lfs=ispc;
if(~is_force_use_no_lfs && strcmp(h.endian,'l'))
    if(should_use_lfs(filename))
        closelfs(fid);
    else
        fclose(fid);
    end
    [fid,msg] = openlfs(filename);
    if (fid < 0)
        error(msg);
    end
    filemode = 'lfs';
end

% if available, get intended stimulus sequence
if ~options.no_intended_times
    stimsOrig=parse_stim_seq(h);
    % TEH: restrict to the subset in which the valve actually changes
    % state (i.e., if the user puts in a "fake transition" which doesn't
    % actually changed the valve #, there will be no record in the
    % stimulus channel.)
    vtindex = find(diff(stimsOrig(:,1)))+1;
    stims = stimsOrig(vtindex,:);
    if(isempty(stims))
        realstim=[];
        return;
    end

    % Convert times to scan numbers
    stims(:,2) = stims(:,2)*h.scanrate;
    wndL=options.wndL * h.scanrate;
    wndR=options.wndR * h.scanrate;
    
    % If in testing mode, resample at correct rate and add to testing_out
    if options.testing_mode
        testing_out.intended_times = round(stims(:,2)./testing_out.resampling_factor);
        testing_out.intended_valvenumbers = stims(:,1);
    end
end

% Determine the stimulus channel
if ~isfield(options,'feedback')
    if(strcmp(key2value(h.wholeheader, 'hardware'), ...
            'dumb-box-0@chanel.wustl.edu'))
        options.feedback = 0;
    else
        options.feedback = 63;
    end
    warning(['warning---assuming ' num2str(options.feedback) ...
        ' is the stimulus channel in timestim().']);
end
feedback_ch_idx=find(h.channels==options.feedback);
if isempty(feedback_ch_idx)
    error(['Feedback channel ' num2str(options.feedback) ' not recorded.']);
end

% Show the progress "bar"
%figProgress=figure;
%tName=['timestim() is processing ' filename ': %d %% done'];
%set(figProgress, 'menubar', 'none');
%tFigPos=get(figProgress, 'position');
%set(figProgress, 'position', [tFigPos(1:2) 500 1]);
%set(figProgress, 'numbertit', 'off')
%set(figProgress, 'name', sprintf(tName, 0));
%shg;

% If in testing mode, go ahead and fill in the resampled voltage trace now,
% since it's too much of a pain to try to integrate the parts already read
% in below (and since it's just testing code, trying to make it clear
% rather than efficient)...
if options.testing_mode
    nScansInMemoryGoal = 1000000;
    nResamplesPerChnk = round(nScansInMemoryGoal/testing_out.resampling_factor);
    chnk_sz = nResamplesPerChnk*testing_out.resampling_factor;;
    nChnk = ceil(h.nscans/chnk_sz);
    resampling_points = (1:nResamplesPerChnk)*testing_out.resampling_factor-(testing_out.resampling_factor-1);
    testing_out.voltage_trace = [];
    for nthChnk = 1:nChnk
        data_start = 1+(nthChnk-1)*chnk_sz;
        data_end = min([(nthChnk*chnk_sz) h.nscans]);
        feedback = readthedata(data_start,data_end,h,fid,filemode,feedback_ch_idx);
        trace = double(feedback);
        if nthChnk == nChnk
            points_to_kill = find(resampling_points>length(trace));
            resampling_points(points_to_kill) = [];
        end
        resampled_trace = trace(resampling_points);
        testing_out.voltage_trace = [testing_out.voltage_trace resampled_trace]; 
    end
end
    
% If get to use intended Times...
if ~options.no_intended_times
    nStims=size(stims, 1);
    realtime=zeros(nStims,1) - 1; % initialized as -1 to ease error detection
    voltages=zeros(nStims,2);      % voltage before,after each transition
    % (at the end we'll do a regression and
    % compare with the expected valve numbers)

    for stimIdx=1:nStims
        % Calculate the expected times
        % (TEH: Note we calculate all of them, because we'll need to look
        % ahead to make sure we don't include part of the next transition)
        if (stimIdx < 2)
            intendedTimes=stims(:, 2);  % Just use the expected values
        elseif (stimIdx == 2)
            intendedTimes = stims(:,2) - intendedTimes(1) + realtime(1);
        else
            % Do a linear regression on transition times, and compute the
            % expected time from that regression
            [slope,offset] = linregress(stims(1:stimIdx-1,2), ...
                realtime(1:stimIdx-1));
            intendedTimes = slope*stims(:,2) + offset;
        end


        % Determine the time window to load
        if (stimIdx == 1)
            data_start=max([0 intendedTimes(stimIdx)-wndL]);
        else
            % TEH: make sure don't include part of previous valve transition
            data_start=max([0 intendedTimes(stimIdx)-wndL realtime(stimIdx-1)+1]);
        end
        if (stimIdx == nStims)
            data_end = min([h.nscans-1 intendedTimes(stimIdx)+wndR]);
        else
            % TEH: don't include next valve transition
            data_end = min([h.nscans intendedTimes(stimIdx)+wndR ...
                intendedTimes(stimIdx+1)-10]);
        end
        if (data_end < data_start)
            error('No available time range for transition');
        end

        % Read the data
        feedback = readthedata(data_start,data_end,h,fid,filemode,feedback_ch_idx);

        % Look for the transition
        [tstep,vfit, ttt, deadloop] = fitstep(double(feedback));

        if(deadloop)
            if options.fix_one_point_causes_of_deadloops_on
                % it seems that in rare cases, one or two single-point blips
                % positioned the right distance away from the actual step
                % can deadloop fitstep.  see if any of these are present
                % using a very dumb check, and, if so, rerun fitstep...
                meanvalue = mean(feedback);
                calcstd = std(feedback);
                ptdev = feedback-meanvalue;
                iwayhigh = find(ptdev>2*calcstd);
                if length(iwayhigh)<3 % not too many points...
                    % and not sequential...
                    ints = iwayhigh(2:end)-iwayhigh(1:(end-1));
                    if ~any(ints<2)
                        for nthwayhigh = 1:length(iwayhigh)
                            if iwayhigh(nthwayhigh)<length(feedback)
                                feedback(iwayhigh(nthwayhigh)) = ...
                                    feedback(iwayhigh(nthwayhigh)+1);
                            else
                                feedback(iwayhigh(nthwayhigh)) = ...
                                    feedback(iwayhigh(nthwayhigh)-1);
                            end
                        end
                    end
                    [tstep,vfit, ttt, deadloop] = fitstep(double(feedback));
                end
                if deadloop % ie, if didn't fix it...
                    error(['dead loop in fitstep() when scan# is: ' num2str(stims(stimIdx,2))]);
                else
                    warning(['removed single-point blips in stimtrace to prevent timestim crashing (scan# approx ' num2str(stims(stimIdx,2)) ')'])
                end
            else
                error(['dead loop in fitstep() when scan# is: ' num2str(stims(stimIdx,2))]);
            end
        end

        % Adjust the position of the step to a real scan number
        realtime(stimIdx) = round(data_start + tstep - 1);
        voltages(stimIdx,:) = vfit;
        % Update the progress meter
        %set(figProgress, 'name', sprintf(tName, round((stimIdx-1)/nStims*100)));
        %drawnow
    end

    % Now we have all the times; let's check to see that the voltages are
    % consistent
    actualVlvs = [[stimsOrig(1,1);stims(1:end-1,1)],stims(:,1)];
    [slope,offset] = linregress(voltages,actualVlvs);
    calcVlvs = round(repmat(slope,size(voltages,1),1).*voltages + ...
        repmat(offset,size(voltages,1),1));
    if any(actualVlvs ~= calcVlvs)
        issame = actualVlvs == calcVlvs;
        issame = issame(:,1).*issame(:,2);
        fractionWrong = 1-sum(issame)/length(issame);
        if fractionWrong > options.calc_vs_intended_tolerance
            error(['Not all valve numbers are consistent: fraction wrong = ' num2str(fractionWrong)]);
        elseif fractionWrong > 0
            warning('There is a discrepancy between calculated and intended valve numbers; however, it is within your specified tolerance.  Please inspect generated figure to ensure you are not missing a hidden recording error....')
            figure
            plot(issame)
            xlabel('Trial Number')
            ylabel('Calculated == Intended?')
            title(['Timestim valvefitting results: Fraction Wrong = ' num2str(fractionWrong)])
        end
    end
    realstim = [[stimsOrig(1,1);stims(:,1);stimsOrig(end,1)],...
        [0;realtime;h.nscans-1]];
    if options.testing_mode
        testing_out.calculated_times = round(realtime./testing_out.resampling_factor);
        testing_out.intended_valvenumbers_matrix = actualVlvs;
        testing_out.measured_voltages_matrix = voltages;
        testing_out.calculated_valvenumbers_matrix = calcVlvs;
        testing_out.voltage2vlvnum_slope = slope(1);
        testing_out.voltage2vlvnum_offset = offset(1);
        testing_out.timestim_final_output = realstim;
    end
else % if not using intended times...
    % Divide into chunks where the max size is a bit smaller than
    % what was given by the user as the min time between valve changes
    min_t_gap = options.min_time_between_valve_changes;
    chnk_sz = .8*min_t_gap*h.scanrate;
    nChnk = ceil(h.nscans/chnk_sz);
    voltages = [];
    change_times = [];
    % Determine the minimum change in voltage reading that should be
    % recognized as a "new" voltage
    if isstr(options.voltage_to_vlvnum_assignments) % If don't yet know the voltages that matter, look for any change that's over .05% of the total spread of voltages used
        % spotcheck a series of chunks of data for the minimum and maximum
        % voltages recorded
        nTo_check = 10;
        time_for_chnk = 10; %in secs
        minVolt = [];
        maxVolt = [];
        start_options = (1:nTo_check)*(h.nscans/(nTo_check+1));
        for n = 1:nTo_check
            data_start = start_options(n);
            data_end = h.scanrate * time_for_chnk + data_start;
            feedback = readthedata(data_start,data_end,h,fid,filemode,feedback_ch_idx);
            trace = double(feedback);
            minVolt = min([minVolt min(trace)]);
            maxVolt = max([maxVolt max(trace)]);
        end
        min_v_gap = (.005) * (maxVolt-minVolt);
    else % If already know the voltages that matter, figure anything below 10% of the smallest change that should matter can be not even looked for
        if size(options.voltage_to_vlvnum_assignments,2) == 2
            voltages = sort(options.voltage_to_vlvnum_assignments(:,2));
            min_v_gap = (.1)*min(voltages(2:end)-voltages(1:(end-1)));
        else
            voltage_mins = sort(options.voltage_to_vlvnum_assignments(:,2));
            voltage_maxs = sort(options.voltage_to_vlvnum_assignments(:,3));
            voltage_gaps = voltage_mins(2:end) - voltage_maxs(1:(end-1));
            min_v_gap = (.1)*min(voltage_gaps);
        end
    end
    % Read in the data in chunks, and look for a transition in each one
    nOverlapScans = 1+options.nToAverage;
    % want to make sure (a) don't just happen to miss
    % miss a valve transition in the space between
    % chunks, and (b) have enough room past a new
    % valve transition to do a bit of averaging if
    % needed
    for nthChnk = 1:nChnk
        % Define the range of data to be read in
        data_start = 1+(nthChnk-1)*chnk_sz;
        data_end = min([(nthChnk*chnk_sz+nOverlapScans) h.nscans]);
        feedback = readthedata(data_start,data_end,h,fid,filemode,feedback_ch_idx);
        trace = double(feedback);
%                 % just for testing
%                 figure(1)
%                 clf
%                 plot(trace)
%                 title(nthChnk)
        % if this is the first chunk, enter the first values...
        if nthChnk == 1
            voltages = mean(trace(1:options.nToAverage));
            change_times = 0;
        end
        % See if the voltage at the end matches the most recent entry into
        % "voltages" (safer than just taking the value of the trace at
        % t=1 because lets you use averaging feature), in which case we move on,
        % or doesn't (by more than the min_v_gap found above), in which case
        % we decide there must be a valve transition here and go about
        % locating it
        v_diff1 = abs(mean(trace((end-options.nToAverage):end))-voltages(end));
        v_diff2 = abs(trace((end-options.nToAverage))-voltages(end));
        % tricky point.  want to use averaging so that we're not
        % fooled by short blips, but we don't want to artificially
        % trigger edge detection by averaging in a large number
        % past our interest zone.  going to bet that our biggest
        % concern is a blip that accidentally just barely triggers
        % an edge detection, and that we're not really at risk for
        % a blip that takes something that otherwise would be an
        % edge and brings it too far down!  so, require BOTH to be
        % "positive" checks
        if (v_diff1 > min_v_gap) & (v_diff2 > min_v_gap)
            % ok, fine, we actually have to do some work now....
            [v,t] = find_v_t(trace,options.nToAverage,voltages(end));
%             % just for testing
%             figure(1)
%             clf
%             plot(trace)
%             hold on
%             plotvertical(t,'r-')
%             title(nthChnk)
            voltages = [voltages v];
            change_times = [change_times (t+data_start)];
            % put checkpoint: make sure alternating signs for voltage
            % changes...
            if length(voltages)>2
                sign1 = (voltages(end)-voltages(end-1)) > 0;
                sign2 = (voltages(end-1)-voltages(end-2)) > 0;
                if sign1 == sign2
                    error('something very wrong!')
                end
            end
        end
    end
    % output histogram of voltages found whether user wants it or not
    fig_handle = figure;
    hist(voltages,1000);
    [ho,xout] = hist(voltages,1000);
    title('Frequency distribution of post-step voltages')
    xlabel('Voltage')
    % Ok, now three possibilities for the actual output...
    if isstr(options.voltage_to_vlvnum_assignments)
        if strcmp(options.voltage_to_vlvnum_assignments,'return_voltages')
            % means output is just supposed to be the vector of voltages alone;
            realstim = [voltages' change_times'];
        elseif strcmp(options.voltage_to_vlvnum_assignments,'guess_mapping')
            % means going to make an educated guess as to what the valve
            % mapping is based on our histogram of voltages.  not terribly
            % efficient, but should only be done once, occasionally, so
            % shouldn't matter...
            % a) Start by determining the distribution of gaps between
            % voltages observed, since this holds most of the
            % information...
            filled_bins = find(ho);
            gaps = filled_bins(2:end)-filled_bins(1:(end-1));
            gaps = gaps(find(gaps>3));
            % b) Now, want to find largest common divisor of all of these
            % gaps.  First see if just the smallest gap itself works
            expected_spread_within_cluster = 5;
            smallest_gap = min(gaps);
            similarly_small_gaps = gaps(find(gaps-smallest_gap<expected_spread_within_cluster));
            smallest_gap_mean = mean(similarly_small_gaps);
            sumError = 0;
            for nthGap = 1:length(gaps)
                sumError = sumError + flexmod(gaps(nthGap),smallest_gap_mean);
                out = [nthGap gaps(nthGap) smallest_gap_mean flexmod(gaps(nthGap),smallest_gap_mean) sumError];%just for testing
            end
            if sumError > expected_spread_within_cluster*length(gaps)
                error('Doesnt look like the first guess of spacing worked, and havent yet implimented the second round; when do, should just progressively divide the smallest gap by 2,3,4 etc until one works')
            else
                slope = 1/(smallest_gap_mean*(xout(2)-xout(1)));
            end
            % c) find offset
            [trash,Imostcommon] = max(ho);
            if Imostcommon > 3
                % this is odd - means, eg, not using first valve for flush
                % solution - or we're way off just in general...
                warning('Something atypical about voltage distribution - please inspect by hand before continuing')
                keyboard
            end
            offset = xout(1);
            % d) turn into valvenumbers!
            calcVlvs = round( (voltages-offset)*slope );
            figure
            hist(calcVlvs,1000)
            xlabel('Valvenumber assigned')
            ylabel('Number of times assigned')
            title('Valvenumber assignments: compare w distribn of post-step voltages')
            realstim = [[calcVlvs 0]' [change_times (h.nscans-1)]'];
            testing_out.timestim_final_output =realstim;
            testing_out.measured_post_step_voltages = voltages;
            testing_out.measured_change_times = change_times;
        end
    else
        % means need to convert output to actual valve numbers
        nVoltageCategories = size(options.voltage_to_vlvnum_assignments,1);
        if size(options.voltage_to_vlvnum_assignments,2) == 2
            % this is going to be inefficient but i don't care right now
            nVoltages = length(voltages);
            setpoints = options.voltages_to_vlvnum_assignments(:,2);
            valveAssignments = [];
            for nthV = 1:nVoltages
                vDiffs = abs(setpoints-voltages(nthV));
                [diff1 iClosest] = min(vDiffs);
                vDiffs(iClosest) = 0;
                [diff2 iNextClosest] = min(vDiffs);
                fractionalDiff = diff1/(diff1+diff2);
                if fractionalDiff>.25
                    error('Valve assignment is ambiguous; please check for errors and/or reconsider voltage tolerances')
                end
                valveAssignments(nthV) = options.voltages_to_vlvnum_assignments(iClosest,1);
            end
            % add lines indicating where the voltage assignment centers
            % are...
            figure(fig_handle)
            hold on
            for nth = 1:nVoltageCategories
                plotvertical(setpoints(nth),'r-')
            end
        else
            nVoltages = length(voltages);
            mins = options.voltages_to_vlvnum_assignments(:,2);
            maxs = options.voltages_to_vlvnum_assignments(:,3);
            valveAssignments = [];
            for nthV = 1:nVoltages
                aboveMin = (voltages(nthV)>mins);
                belowMax = (voltages(nthV)<maxs);
                iRight = find(aboveMin & belowMax);
                nRight = length(iRight);
                if nRight ~= 1
                    error('Either no valve assignment or more than one valve assignment...')
                end
                valveAssignments(nthV) = options.voltages_to_vlvnum_assignments(iRight,1);
            end
            % add indicators of how the voltages were sorted to valvelabels
            figure(fig_handle)
            hold on
            for nth = 1:nVoltageCategories
                middle = mean([mins(nth) maxs(nth)]);
                plotvertical(middle,'r-')
                plotvertical(mins(nth),'g:')
                plotvertical(maxs(nth),'g:')
            end
        end
        realstim = [valveAssignments' change_times'];
    end
end

if options.testing_mode
    save timestim_testing_out testing_out
end







%%--------------------------------------------------------

function [feedback] = readthedata(data_start,data_end,h,fid,filemode,feedback_ch_idx)

% Read the data
if strcmp(filemode,'lfs')
    tAllChData = readint16lfs(fid,h.numch,[data_start data_end], ...
        h.headersize);
else
    fseek(fid,h.numch*data_start*2+h.headersize,'bof');
    tAllChData = fread(fid,[h.numch,data_end-data_start+1],'int16');
end
feedback=tAllChData(feedback_ch_idx, :);


%----------------------------------------------------------

function [v,t] = find_v_t(trace,nToAverage,v1)

t1 = 1;
t2 = length(trace)-nToAverage;
v2 = trace(t2);
time_gap = t2-t1;
while time_gap > nToAverage
    tmid = round(time_gap/2)+t1;
    vmid_plus = mean(trace(tmid:(tmid+nToAverage)));
    start_vmid_minus = max([(tmid-nToAverage) (.6)*(nToAverage)]); % otherwise if transition very close to start could end up trying to index times before start of this chunk
    vmid_minus = mean(trace(start_vmid_minus:tmid));
    dva = abs(vmid_plus-v1);
    dvb = abs(v2-vmid_minus);
    if dva>dvb
        t2 = tmid;
        v2 = vmid_plus;
    elseif dva < dvb
        t1 = tmid;
        v1 = vmid_minus;
    else
        error('Hard to see how this could happen unless the interval shouldnt have been picked to begin with')
    end
    time_gap = t2-t1;
end
t = t2;
v = v2;



