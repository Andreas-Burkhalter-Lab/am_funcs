% Functions for analyzing and plotting analytical chemistry data
%
% CENTROID-MODE DATA:
% Loading centroid mass spectrum data:
%   msload                     - Load mass spec data and parse into peaks.
%   ms2peaks                   - Parse mass spec data into peaks.
%   mschoosetrials             - Load mass spec data, let user choose trials (GUI).
%   isotopes                   - load NIST information about masses and abundance of isotopes
%   msbrowse                   - interactively inspect a mass-spec file
%   msload_mzXML               - load mass spectrometry data from mzXML format
%   msload_netcdf              - load mass spec .cdf files from Gross lab
%   ms_set_gradient_from_table -save gradient state in each scan
%   resample_gradient          -set resampling times to yield a common linear gradient
%   scan2matrix                - convert mass spec data to sparse matrix format

% Analyzing centroid mass spectrum data:
%   mschoosepeak    - Pick biggest unassigned peak.
%   msmatchpeak     - Find a given peak in all fractions.
%   mslabelpeaks    - Label a particular peak in all fractions.
%   msgetcompounds  - Pick out the largest peaks in all fractions.
%   formulazmass    - calculate mass of the "base" compound with given formula
%   formula_chisq   - calculate goodness-of-fit for candidate formulas to isotopologue data
%   isotopologues   - generate the complete list of isotopologues for given molecular formula
%   mass2formula    - generate list of possible molecular formulas producing a given mass
%   msdomain_stats.m
%   mscollect_points - options fields: mz_max, mz_min
%   msdecimate       - ptions fields: n_bins_mZ
%   mspeaks_by_imflow.m 
%   mz2indx          -convert m/z values to an index
%   peaks2compounds  -convert profile-mode data to centroid-compatible
%   string2formula   -calculate a numeric "formula" given a chemical formula string



% PROFILE-MODE DATA:

%   msprof_bgfg                            - load mass spec profile data, with background subtraction
%   msprof_bgfgplot.m
%   msprof_charge_isotopes_linear          - determine charge and element composition from isotopologues
%   msprof_compounds_pairwise_polt         - compare concentrations of compounds between samples
%   msprof_chooase_peak_timecourses_gui    -for msprof_choose_peak_timecourses_gui.fig
%   msprof_concentrations                  -measure compound concentrations from NNMF and timeSpans
%   msprof_consolidate_nnmf                - convert structure to matrix format
%   msprof_consolidate_nnmf.m              - convert structure to matrix format
%   msprof_correlate_dr                    -correlate physiological response with concentration
%   msprof_correlate_gui                   -correlate neuronal responses to ligand abundance
%   msprof_correlate_gui_temp.m            -correlate neuronal responses to ligand abundance
%   msprof_extract_compounds_set           - process a set of mass spec profiles
%   msprof_isotope_profile                 -create matrices facilitating isotopologue fitting
%   msprof_load_set                        -load a set of mass spec experiments to matrix format
%   msprof_monotonic_dr                    -compare concentration & firing rate data using a monotonic model
%   msprof_nnmf                            -decompose a set of mass spec runs into m/z peaks
%   msprof_nnmf_pairwise_plot.m 
%   msprof_plot_isotopolgoues              -compare measured and expected abundances of isotopologues
%   msprof_run                             - Load & extract peaks from profile mass spec data
%   msprof_temporal_register_gui           -correct errors in gradient timing


% Plotting mass spectrum data:
%   msplotcompound  - Plot concentration vs. fraction # + GUI m/z spectrum.
%   formula_gui     - inspect mass spec data related to chemical formula
%   msplot.m

% Correlating to firing rates:
%   mscorrelatedr   - Compare concentrations & firing rates.
%   
%
% Helper functions
%   mschoosetrials_gui   - GUI for a single fraction.
%   msprof_parse         - Load profile-mode mass spec data
%   msprof_factor        - Separate into background + signal
%   msprof_findpeaks     - Identify peaks 
%   msprof_peakamp       - Quantify the size of each peak













