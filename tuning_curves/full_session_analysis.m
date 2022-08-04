 %%% perform rf mapping and tuning curves analysis for each roi
 %       pars = pars struct for get_tuning_curves
 % updated 2019-1-25 on thermaltake
 
 function tuningdat_out = full_session_analysis(input_files,pars)
 
 %% analyze roi responses
 % list all sub-sessions in the concatenated run; will be used to eliminate extra scope pulses in sessions that used multiple planes
 pars.all_scopetiming_files = {input_files.scopetiming_file_rf, input_files.scopetiming_file_dark, input_files.scopetiming_file_tuning};
% tuning curves 
tuningpars = get_tuning_curves(input_files.dff_file, input_files.triggertiming_file, input_files.stimdata_file_tuning, input_files.scopetiming_file_tuning, input_files.regops_file,...
    input_files.pupil_data_file, input_files.patchdata_file, input_files.patch_reg_file, pars);

%% combine rf analysis and tuning curves analysis
if isempty(input_files.scopetiming_file_rf) %%% if we're only analyzing sf/tf/orient/size tuning, and not rf mapping
    tuningdat_out = tuningpars; 
elseif ~isempty(input_files.scopetiming_file_rf)
    rfdat = get_tuning_curves(input_files.dff_file, input_files.triggertiming_file, input_files.stimdata_file_rf, input_files.scopetiming_file_rf, input_files.regops_file,...
        input_files.pupil_data_file, input_files.patchdata_file, input_files.patch_reg_file, pars);
    rf_and_stim_overlap() % find out the overlap of the tuning stim and each roi's receptive field
    tuningdat_out = struct;
    tuningdat_out.input_files = input_files;
    % variables to combine from rf struct and tuning struct
    variables_to_combine = {'meanImage_pre_rotate', 'regops_file', 'nplanes_in_session_total','iplane','triggertiming_file','dff_file',...
        'patchdata_file','patchdata','registration_file','reg_struct'};
    for i = 1:length(variables_to_combine)
        thisvarname = variables_to_combine{i};
        try ismatch = all(tuningpars.(thisvarname) == rfdat.(thisvarname));
            if ~ismatch
                error(['tuning and rf structs don''t match in field ''' thisvarname ''])
            end
        end
        if isfield(tuningpars,thisvarname) % if this variable was saved into tuningpars
            tuningdat_out.(thisvarname) = tuningpars.(thisvarname); 
            tuningpars = rmfield(tuningpars, thisvarname);
            if ~isempty(input_files.scopetiming_file_rf)
                rfdat = rmfield(rfdat, thisvarname);
            end
        end
    end
    tuningdat_out.stimpars = tuningpars.stimpars; %%% add to top struct for compatibility with various funcs
    % combine roi tables... get rid of variables that cause errors for join function
    rfdat.tuning.roi_image_prereg = [];
    rfdat.tuning.roi_image_pre_rotate = [];
    rfdat.tuning.centeryx_pre_rotate = [];
    if  any(ismember(rfdat.tuning.Properties.VariableNames,'roi_image_reg'))
        rfdat.tuning.roi_image_reg = [];
    end

    tuningdat_out.tuning = join(tuningpars.tuning, rfdat.tuning);
    % combine scope_events tables
    tuningdat_out.scope_events = tuningpars.scope_events; 4
    tuningdat_out.scope_events = [tuningdat_out.scope_events, rfdat.scope_events(:,{'stim_row','stim_column'})]; % add rf stim data to table
    tuningdat_out.scope_events = movevars(tuningdat_out.scope_events,{'stim_row','stim_column'},'After','stim_diam_deg');
    tuningdat_out.scope_events.stim_present = tuningpars.scope_events.stim_present | rfdat.scope_events.stim_present; % if either stim type is present, note it

    % clean up extra fields
    if all(all(isnan(tuningdat_out.tuning.centeryx_pre_rotate))) % if we didn't rotate image
        tuningdat_out.tuning.centeryx_pre_rotate = [];
        tuningdat_out.tuning.roi_image_pre_rotate = [];
    end
    if any(ismember(tuningdat_out.tuning.Properties.VariableNames,'nearest_patch_center_yx'))  &&  all(isnan(tuningdat_out.tuning.nearest_patch_center_yx))
        tuningdat_out.tuning.nearest_patch_center_yx = [];
    end
    tuningpars = rmfield(tuningpars,{'tuning','scope_events'});
    rfdat = rmfield(rfdat,{'tuning','scope_events'});
    tuningdat_out.tuningpars = tuningpars;
    tuningdat_out.rfdat = rfdat; 
end

%% locomotion correlation
if ~isempty(input_files.scopetiming_file_dark)
    tuningdat_out = locomotion_analysis(tuningdat_out); 
end
