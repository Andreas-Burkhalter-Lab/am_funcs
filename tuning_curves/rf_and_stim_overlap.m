%%% find the proportion of the tuning stimulus which overlaps with each ROI's receptive field
% called by full_session_analysis
%
%%% updated 2020/04/25 on thermaltake

if isfield(tuningpars.stimpars,'stim_area') && pars.fit_tuning_functions % if the tuning stimulus area was saved and rf image was computed
    stim_area = tuningpars.stimpars.stim_area; 
    npix_stim_area = nnz(stim_area);
    nrois = height(rfdat.tuning); 
    for iroi = 1:nrois
        this_roi_rf_image = rfdat.tuning.rf_fit{iroi}.rf_image;
        npix_overlap = nnz(this_roi_rf_image & stim_area); % number of screen pix shared by rf and stimulus
        rfdat.tuning.stim_on_rf(iroi) = any(npix_overlap); % true if the stimulus covered any of this roi's RF
        rfdat.tuning.stim_overlap_frac(iroi) = npix_overlap / npix_stim_area; % proportion of stim that was in this roi's RF
        % get proportion of rf that was covered by the stim; approximate because it doesn't account for parts of the rf that are off the stimulus screen
        rfdat.tuning.rf_overlap_frac_aprx(iroi) = npix_overlap / nnz(this_roi_rf_image); 
    end
end
