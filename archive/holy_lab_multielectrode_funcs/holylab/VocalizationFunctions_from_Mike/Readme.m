% This is a list of all the scripts which are currently in use.



% bdurcorr: 	Looks for within-bout correlations 1,2,3 neighbors away of whistle duration. Requires a minimum of 5 data points for each comparison (Meaning 5 x1vs.x2, x1 vs. x3, x1 vs. x4).
% 
% bootsngsci:	Outputs means and bootstrapped 95% confidence intervals for a cell array x (1x257 power values), and N bootstrap population size.
%
% boutdetector:	Finds bout boundaries given cutoff determined from pause distribution.
%
% bpausecorr:	Looks for within-bout correlations 1,2,3 neighbors away with a minimum of 5 data points for each comparison of intrabout pause lengths.
%
% cilines:		Edited per use, draws up-and-down lines for 95% cis determined by bootsngsci.
%
% fakesng:		Not fully implemented, but makes an artificial sonogram of a single whistle type based on a model whissnip for that whistle type, with randomly distribution interbout and intrabout lengths to build up a 3 minute SNG. The interbout and intrabout lengths are determined by random length from within the ranges determined for the Celf6 2012 data. There was no attempt to respect correlations of one whistle to the next.
%
% logpowerspectrum:	Calculates 1x257 with log10 power values, rather than relative (fractional) power values.
%
% make_whistles: Runs whissnip, whistimes, whisparams for all SNGs in a directory and saves output as *.mat files.
%
% makea3dhist: Makes 3D histogram with 10,000 bins (100x100) for data of <t,t+1msec> whistle peak frequencies.
%
% manualthresholds: Environment for plotting 20 sec samples from each SNG, plotting max power in each msec time bin over all bins for selecting background threshold. Requires an SNG built with no threshold (sngparms.threshold = 0;)
%
% params_explained: Comment file with explanation of each structure field in the output of whisparams
%
% powerspectrum: Calculates relative (fraction of total) power in each 1x257 frequencies for each SNG.
%
% pupsong_pipeline: Not intended to be run all at once, the full pipeline for pupsong analysis as of 08/19/2013. Calls other scripts in this Readme file, but there are some modules which have yet to be edited as their own functions. Run each part separately, don't run as a single command! This was built on the summer 2013 178C pup data, but now is a template. Save as a separate file under the directory of interest and then run the modules, in case you have to make edits.
%
% Readme: This readme file.
%
% scatcat: Plots data as dots scattered randomly with a certain 'factor' of scatter around an x value.
%
% semgraph: Plots a line and bars for mean ± standard error on any graph.
%
% shuffle: Shuffles a 1xn vector x N times.
%


