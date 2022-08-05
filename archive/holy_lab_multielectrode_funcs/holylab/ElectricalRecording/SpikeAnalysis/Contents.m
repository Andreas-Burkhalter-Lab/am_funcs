 Routines for the analysis of electrophysiological recordings
%
% ephys structure utilities:
%  copyotherfields.m - copy some fields from one structure to a new one
%  ephys          - Describes the fields of the ephys structure.
%  ephysfromai    - Set up an ephys structure from a AI headers.
%  ephysjoincells - Join two or more cells together.
%  ephyssubchan   - Select a subset of channels.
%  ephyssubcell   - Select a subset of cells.
%  ephyssubrange  - Select a set of time intervals.
%  ephysfetch     - Read ephys data from files.
%  ephystofetch   - Determine which data need fetching.
%  ephystag       - Name epochs according to strings, durations, etc.
%  ephys_cass_fetch.m - a temporary means for importing CASS sorting results
%  ephyschannelbrowser.m - examine responses in single file
%  ephysconvert.m - convert older data formats into ephys formats
%  ephysfilebrowser.m - GUI for browsing files and channels
%  ephysplottimecourse.m - show PSTH for a whole experiment
%  ephysrethresh.m - apply a new (higher) threshold for spike detection
%  ephysvalidatetimerange.m - check for consistent time range across repeats


%
% Analysis functions:
%  ephysbinspikes - Bin spikes in equal time windows.
%  ephyspsth      - Calculate the PSTH across repeats.
%  deltarate      - Compute firing rate change upon stimulation.
% ddeltarmaxbytime.m - compute rate difference upon stimulation with time onset consideration
% deltarmax_calculation.m - computing maximum firing rate in a response
% deltarmax_HAAmod.m - computing maximum firing rate in a response, modified by HAA
% deltarmax.m- compute maximum rate difference upon stimulation
% deltar_statistics.m - compute mean and sem of change in firing rate
% drmax_acrosstrials.m - compute rate difference upon stimulation. Finds max dr across all trials
% findst_time.m - this function will pool the spike times across the array to figure out
% rmax_seriesHAAmod.m - compute neural response to concentration series
% rmax_series.m - compute neural response to concentration series

%
% Simulated cells:
%  ephyssimulatedspikes - simulated spike trains using a firing rate model.
%
% Plotting functions:
%  drgui.m 	   - a graphical browser for electrophysiology results
%  ephysplotparams - Describes the fields which control the plot.
%  ephysplot       - Plot data in an ephys structure.
%  ephysgui        - Plot ephys data and set up GUI.
%  deltarategui    - Plot firing rate changes upon stimulation.
%  stim_spike      - Analyze spike waveforms & stimulus responses w/o sorting.
% mdexplore_ephys.m - connect an mdexplore plot to ephys data
%
% Variable time-window binning:
%  SplitEvenly       - Bins with equal numbers of spikes/bin.
%  SplitGaussEven    - Equal #/bin, # chosen by skewness threshold.
%  SplitGauss        - Choose bin boundaries to minimize skewness.
%
% Modelling:
%  conccompare     - Fit spike data to a ligand binding model.
%  conccompareplot - Graphical output from the above analysis.
%  concratio       - Summarize effective concentration ratios for many cells.
%  concchisq       - Inner loop for the fitting process in conccompare.
%  fit1storderss.m - Fit steady-state responses to Michaelis-Menten function
%  rate1storder.m  - Compute firing rates from a 1st order binding model
%
%
% Interval utilities:
%  intervalsfromstim - Choose time intervals based on stimulus.
%  timesininterval   - Select those times falling in specified intervals.
%  ephysvalidatetimerange - Enforce consistent time range across repeats.
%  organize_by_cycle - Collect presentations into single-repeat cycles.
%  refractory_violation_stats.m - statistics on refractory violations
%
% Julian Meeks' code:
%    calc_delta_r_from_ephys.m - Calculates firing rate changes from ephys data
%    calc_fem_selectivity.m - Calculates selectivity from drF and drM
%    calc_fem_selectivity_stderr.m - Calculates std of selectivity from drF and drM
%    calc_lifetime_sparseness.m - Calculates lifetime sparseness of responses
%    calc_lifetime_sparseness_neg.m - Calculates lifetime sparseness of negative responses
%    calc_lifetime_sparseness_pos.m - Calculates lifetime sparseness of positive responses
%    calc_mix_subadd.m - Calculates subadditivity from drF, drM, drMix and base_r
%    calc_mix_supp.m - Calculates mixture suppression from drP, drMix and base_r
%    calc_selectivity_from_ephys.m - calc_selectivity_from_ephys (JPM, script only at this time)
%    meea_apply.m 	- Applies the most recent version of the seea_analysis to the directories in directory_list
%    meea_baserate 	- Calculates baseline firing rates for each cell in input
%    meea_deltar.m 	- calculates delta_r values over the designated time range(s)
%    meea_maxsum_timeweight.m - modifies deltar_max vals based on early/late tmax
%    meea.m 		- MEAA is a shorthand for MULTI_ELECTRODE_EPHYS_ANALYZE: analyzes current directory with default values
%    meea_maxsum 	- Chooses finds tmax and rmax for multielectrode ephys analysis
%    meea_maxsum_tstart 	- Sets time offsets for multielectrode ephys analysis
%    meea_maxsum_valvegroup.m - Parses valvelabels to group compounds for maxsum_monotonic
%    meea_resample 	- Adds a resampled firing rate at the desired rate(s)
%    multi_electrode_ephys_analyze.m - Analyzes multi electrode data from a directory
%    strip_se.m 		- Parses s%e% notation to return the %values
%    jpm_ephys_script1.m - Julian's 11-29-07 script to extract and plot the things for the NRSA Dec 2007 deadline
%    jpm_ephys_script2.m - Julian's 11-30-07 script to extract and plot the things for the NRSA Dec 2007 deadline
%    jpm_ephys_script3.m - Julian's 2008_01_24 script to extract and plot a running tab of delta_r values from the AOB sulfated steroid project
%    jpm_ephys_script4.m - Julian's 2008_04_29 precursor to a more robust analysis package for ephys data
%    seea_apply.m - applies the most recent version of the seea_analysis to the directories in directory_list
%    seea_deltar.m - calculates delta_r values over the default time range
%    seea_discrim.m - calculates d' (discriminability) index from ephys or analysis struct
%    seea_lifetime_sparseness.m - calculates lifetime sparseness for sulfated steroids in a given experiment
%    seea.m - SEAA is a shorthand front for SINGLE_ELECTRODE_EPHYS_ANALYZE
%    seea_mix.m - calculates differences between mixtures of identical stimuli
%    seea_plot_stats.m - creates a 3-d plot of p values in an analysis struct
%    seea_selectivity.m - calculates gender selectivity on all concentrations (dilutions) of male and female urine
%    single_electrode_ephys_analyze.m - analyzes single electrode data from directory
%    test_seea_func.m - Test function
%
% Outdated:
%  BinSpikes  - Bin spikes using old cell array format.
%  PSTH       - Compute PSTH using cell array format
%
