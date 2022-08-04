%%% load tuningdata from multiple planes and concatenate all rois in a single table
%
% [roitable, filetable] = analyze_multi_planes(filetable_in, excel_rows)
%
% filetable_in = master excel file listing tuning and patchdata files to analyze, or table variable
%%%%% updated 2021-1-15

function [roitable, filetable] = analyze_multi_planes(filetable_in, excel_rows)

rerun_roi_patch_location = 0; %%% rerun get_roi_patch_location.m on each tuningdat struct
rerun_pakan_locmod = 0; %%% rerun pakan_locomotion_modulation.m on each tuningdat struct
get_loc_vs_rest_tuning = 0; %% get tuning data, including HWHM, for high locomotion vs. low locomotion trials
get_discrimination_index = 1; %% get DI for sf, tf, orient (Gao et al. 2010 'Parallel input channels...')

% stats for computing RF locations
screenwidth_pix = 1920; 
screenheight_pix = 1080; 

if ~istable(filetable_in) % if filetable_in is a filename rather than an already-processed table variable
    filetable = excel2table(filetable_in); % load master file list
    excel_rows = vardefault('excel_rows',2:height(filetable)+1);
    matlab_table_rows = excel_rows - 1; % compensate for variable names row taking up one place in excel file
    filetable = filetable(matlab_table_rows,:); % only keep rows specified by excel_rows
    filetable.analyze_plane(isnan(filetable.analyze_plane)) = 0; 
    filetable = filetable(logical(filetable.analyze_plane),:); % keep only planes that have been marked for analysis
    filetable.analyze_plane = []; 
elseif istable(filetable_in)
    filetable = filetable_in;
end
    
nfiles = height(filetable);

    
%%% get roi patch location data, concatenate in single table
wbar = waitbar(0,'Gathering tuning and patch location data...');
for ifile = 1:nfiles
    tuning_file = [filetable.directory{ifile}, filesep, filetable.tuning_file{ifile}];
    patch_file = [filetable.patch_file{ifile}];
    reg_to_patches_file = [filetable.reg_to_patches_file{ifile}];
    tuningdat = load(tuning_file,'tuningdat'); tuningdat = tuningdat.tuningdat;
    nrois = height(tuningdat.tuning); 
    pars_to_analyze = {'sf', 'tf', 'orient'}; % variables to analyze for DI and loc_vs_rest
    
    % re-run cell alignment with M2 quantiles
    if rerun_roi_patch_location
        tuningdat = get_roi_patch_location(tuningdat, patch_file, reg_to_patches_file);
        delete([tuning_file, '.mat'])
        save(tuning_file,'tuningdat')
    end
    
    % re-compute locomotion modulation index
    if rerun_pakan_locmod
        tuningdat.tuning = pakan_locomotion_modulation(tuningdat.tuning, tuningdat.scope_events, tuningdat.tuningpars.scopetiming); 
        delete([tuning_file, '.mat']);         
        save(tuning_file,'tuningdat')
    end
    
    
    %%%% get discrimination indices (Gao et al. 2010, Ji et al. 2015)
    %%%% [Gao et al. equation has typo, correct version is in Ji et al. 2015]
    % DI = Rmax ? Rmin/((Rmax ? Rmin) + 2 ?SSE/(N ? M))
    % Rmax is the mean response to the most effective stimulus
    % Rmin is the response to the least effective stimulus
    % SSE denotes the sum squared error around the mean responses
    % N is the number of observations (trials)
    % M is the number of stimulus values tested
    if get_discrimination_index
        for ipar = 1:length(pars_to_analyze)
            thispar = pars_to_analyze{ipar};
            tuningdat.tuning.([thispar, '_di']) = NaN(height(tuningdat.tuning),1);
            for iroi = 1:nrois % for every cell from this session
                resptable = tuningdat.tuning{iroi, [thispar, '_trials']}{1};
                best_row = find(tuningdat.tuning{iroi, [thispar, '_best_tested_val']} == tuningdat.tuning{iroi, [thispar, '_trials']}{1}{:,[thispar]});
                Rmax = resptable.resp_mean(best_row); 
                [Rmin, worst_row] = min(resptable.resp_mean); % least effective stimulus                
                Ntrials = nnz(~isnan(resptable.resp(:))); % total analyzable trials
                Mvals = height(resptable); % number of stim parameter values test
                 % sum square error computation... adapted from Compute_DDI.m by Greg DeAngelis [used in Ji et al. 2015]
                for istimval = 1:Mvals % get sum squared error for each stimulus value for this ROI
                    resptable.sse(istimval) = sum( [resptable.resp(istimval,:) - resptable.resp_mean(istimval)].^2 );
                end
                SSE = sum(resptable.sse);
                tuningdat.tuning{iroi, [thispar, '_di']} = [Rmax-Rmin] / ([Rmax-Rmin] + 2*sqrt(SSE)/[Ntrials-Mvals]); % Ji et al. 2015
                clear resptable best_row Rmax Rmin worst_row Ntrials Mvals SSE
            end
        end
        delete([tuning_file, '.mat']);         
        save(tuning_file,'tuningdat') % save computed DI results into the subject's tuningdat file
    end
    
    % get tuning data, including HWHM, for high locomotion vs. low locomotion trials
    if get_loc_vs_rest_tuning 
        plotops.show_plot = 0; 
        plotops.get_hwhm = 1; 
        for ipar = 1:length(pars_to_analyze)
            thispar = pars_to_analyze{ipar};
            tuningdat.tuning.([thispar, '_fitresults_low_locm']) = cell(height(tuningdat.tuning),1);
            tuningdat.tuning.([thispar, '_fitresults_high_locm']) = cell(height(tuningdat.tuning),1);
        end
        for iroi = 1:nrois % for every cell from this session
            for ipar = 1:length(pars_to_analyze)
                thispar = pars_to_analyze{ipar}; 
                fitresults = plot_tuning_curve(tuningdat, iroi, thispar, [], plotops);
                tuningdat.tuning{iroi, [thispar, '_fitresults_low_locm']} = {fitresults.fitresults_low_locm};
                tuningdat.tuning{iroi, [thispar, '_hwhm_low_locm']} = fitresults.fitresults_low_locm.hwhm;
                tuningdat.tuning{iroi, [thispar, '_fitresults_high_locm']} = {fitresults.fitresults_high_locm};
                if ~isfield(fitresults.fitresults_high_locm, 'hwhm') || isempty(fitresults.fitresults_high_locm.hwhm)
                    tuningdat.tuning{iroi, [thispar, '_hwhm_high_locm']} = NaN; % if there weren't enough high loc trials to compute hwhm
                else % if there were enough high loc trials to compute hwhm
                    tuningdat.tuning{iroi, [thispar, '_hwhm_high_locm']} = fitresults.fitresults_high_locm.hwhm;
                end
            end
        end
        delete([tuning_file, '.mat']);         
        save(tuning_file,'tuningdat')
    end
    tuningdat.tuning = [repmat(filetable(ifile,{'sub','day','plane','layer'}),height(tuningdat.tuning),1), tuningdat.tuning]; % add session data to roi table
    
    % account for missing responsivity variables
    if ifile==1
        roitable = tuningdat.tuning;        
    else
        if any(ismember(roitable.Properties.VariableNames,'rspv_pval')) && ~any(ismember(tuningdat.tuning.Properties.VariableNames,'rspv_pval'))
            tuningdat.tuning.rspv_pval = NaN(height(tuningdat.tuning),1);
        end
        if any(ismember(roitable.Properties.VariableNames,'rspv')) && ~any(ismember(tuningdat.tuning.Properties.VariableNames,'rspv'))
            roitable.rspv = [];
        end
        if ~any(ismember(roitable.Properties.VariableNames,'rspv')) && any(ismember(tuningdat.tuning.Properties.VariableNames,'rspv'))
            tuningdat.tuning.rspv = [];
        end
        roitable = [roitable; tuningdat.tuning];
    end
    try waitbar(ifile/nfiles,wbar); end
    clear tuning_file patch_file reg_to_patches_file thisdat tuningstruct
end
roitable.rf_size_yx_deg(roitable.rf_size_yx_deg==0) = NaN; %% make sure that rf size was not falsely assigned as zero

roitable.rf_center_pix_y = roitable.rf_center_yx_pix(:,1); % separate rf center vars to make easier to analyze
roitable.rf_center_pix_x = roitable.rf_center_yx_pix(:,2); % separate rf center vars to make easier to analyze
roitable.rf_offcenter_pix_y = abs(roitable.rf_center_yx_pix(:,1) - screenheight_pix/2); % distance in pix of rf from screen center
roitable.rf_offcenter_pix_x = abs(roitable.rf_center_yx_pix(:,2) - screenwidth_pix/2); % distance in pix of rf from screen center
try close(wbar); end





