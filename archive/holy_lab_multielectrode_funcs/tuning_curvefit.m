%TUNING_CURVEFIT: Fit curve to tuning data for parameter of interest.
%%% Called by tuning_main.m
% updated 1/8/16

%% Get tuning data
% If spikes for a given trial are missing, leave as NaN.
% Make table for tuning data
tuningParVals = getfield(stimpars,[stimpars.tuningParameter '_vals']); 
nTuningParVals = length(tuningParVals); % number of values the stimulus parameter of interest takes on
tuning = dataset(tuningParVals','VarNames',{'parVal'});
tuning.meanspikes = NaN(nTuningParVals,1);
tuning.spikes = NaN(nTuningParVals,stimpars.repetitions);
tuningParValSequence = double(stimrec(:,stimpars.tuningParameter)); % tuning parameter sequence from trialdata

% Organize spike responses by tuning parameter value.
for valueInd = 1:nTuningParVals
    thisvalue = tuningParVals(valueInd);
    match = tuningParValSequence == thisValue; 
    tuning.spikes(valueInd,:) = [double(trialdata.spikes(match))]; % all spike responses for this param value
    tuning.meanspikes(valueInd) = nanmean(tuning.spikes(valueInd,:)); % mean spike responses for this param value
end
    

%% Fit curve to tuning data
parValsToPlot = linspace(min(tuningParVals), max(tuningParVals), n_param_vals_forplotting); 
switch stimpars_out.tuningParameter
    case 'orient' % fit von mises... should be sum of 2 von mises (see Gao)
        [orient_fitPeak orient_fitKappa] = circ_vmpar(stimpars.orient_vals,double(tuning.meanspikes));
        [orient_fitSpikesToPlot orient_valsToPlot] = circ_vmpdf(parValsToPlot, orient_fitPeak, orient_fitKappa); 
        tuningVarsToSave = {'orient_fitPeak','orient_fitKappa','orient_fitSpikesToPlot','orient_valsToPlot'};
    case 'sf' % fit Gaussian - get wq/Gao function to fit/
        [sf_fitSigma sf_fitPeak sf_spikesAtFitPeak] = mygaussfit(stimpars.sf_vals,double(tuning.meanspikes));
        sf_valsToPlot = parValsToPlot;
        sf_fitSpikesToPlot = sf_spikesAtFitPeak * exp( -(sf_valsToPlot-sf_fitPeak).^2 / (2*sf_fitSigma^2) );
        tuningVarsToSave = {'sf_fitSigma','sf_fitPeak','sf_spikesAtFitPeak','sf_valsToPlot','sf_fitSpikesToPlot'};
    case 'tf' % fit Gaussian - get wq/Gao function to fit/
        [tf_fitSigma tf_fitPeak tf_spikesAtFitPeak] = mygaussfit(stimpars.tf_vals,double(tuning.meanspikes));
        tf_valsToPlot = parValsToPlot;
        tf_fitSpikesToPlot = tf_spikesAtFitPeak * exp( -(tf_valsToPlot-tf_fitPeak).^2 / (2*tf_fitSigma^2) );
        tuningVarsToSave = {'tf_fitSigma','tf_fitPeak','tf_spikesAtFitPeak','tf_valsToPlot','tf_fitSpikesToPlot'};
    case 'diam' % fit Gaussian - get wq/Gao function to fit/
        [diam_fitSigma diam_fitPeak diam_spikesAtFitPeak] = mygaussfit(stimpars.diam_vals,double(tuning.meanspikes));
        diam_valsToPlot = parValsToPlot;
        diam_fitSpikesToPlot = diam_spikesAtFitPeak * exp( -(diam_valsToPlot-diam_fitPeak).^2 / (2*diam_fitSigma^2) );
        tuningVarsToSave = {'diam_fitSigma','diam_fitPeak','diam_spikesAtFitPeak','diam_valsToPlot','diam_fitSpikesToPlot'};
end

%% Save results
eval(['tuning_' stimpars.tuningParameter '=tuning']); % rename with parameter of interest
clear tuning

