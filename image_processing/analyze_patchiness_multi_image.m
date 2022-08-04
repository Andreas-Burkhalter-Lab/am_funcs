%%% analyze_patchiness_multi_image
%       results = analyze_patchiness_multi_image(excel_filelist_name)
%
%%% input the name of an excel spreadsheet containg a filelist with headers: 
%       m2file,	comparison_file, zoom, scope, roifile,	baselineFile
%
%
% updated 2020/5/11

function results = analyze_patchiness_multi_image(excel_filelist_name)

results = readtable(excel_filelist_name);
nimages = height(results);
results.analysis = cell(nimages,1); 

% perform patchiness analysis for all image sets, store restuls in filelist
for jimage = 1:nimages
    if all(isnan(results.roifile)) || isempty(results.roifile{jimage})
        this_roi_file = [];
    else
        this_roi_file = results.roifile{jimage};
    end
    if all(isnan(results.baselineFile)) || isempty(results.baselineFile{jimage})
        this_baseline_file = [];
    else
        this_baseline_file = results.baselineFile{jimage};
    end
    results.analysis{jimage} = analyze_patchiness(... % run analysis on this image set
        results.m2file{jimage}, results.comparison_file{jimage}, this_roi_file,...
        this_baseline_file, results.zoom(jimage), results.scope{jimage});
    results.patchInterpatchRatio(jimage) = results.analysis{jimage}.permutation_test.patchInterpatchRatio;
    results.shufflemean_patchInterpatchRatio(jimage) = results.analysis{jimage}.permutation_test.shufflemean_patchInterpatchRatio;
    results.pval_permutation_test(jimage) = results.analysis{jimage}.permutation_test.pval_permutation_test;
    results.shuffledist_patchInterpatchRatio{jimage} = results.analysis{jimage}.permutation_test.shuffledist_patchInterpatchRatio;
end

% combine results across images
results.m2file{nimages+1} = 'combined'; results.comparison_file{end} = 'combined'; results.scope{end} = 'combined';
results.patchInterpatchRatio(end) = mean(results.patchInterpatchRatio(1:end-1)); % average across image sets
results.shufflemean_patchInterpatchRatio(end) = mean(results.shufflemean_patchInterpatchRatio(1:end-1));  % average across image sets
all_shuffledist_patchInterpatchRatio = [];
for jimage = 1:nimages % create a distribution which is the mean of distributions from all image sets
    all_shuffledist_patchInterpatchRatio = [all_shuffledist_patchInterpatchRatio; [results.shuffledist_patchInterpatchRatio{jimage}]'];
end    
results.shuffledist_patchInterpatchRatio{end} = mean(all_shuffledist_patchInterpatchRatio,1); % create 'grand average' shuffled distribution
results.pval_permutation_test(end) =... %%% get pval for combined p:ip ratio differing from 1
    length(find(abs(results.shuffledist_patchInterpatchRatio{end}-1)>abs(results.patchInterpatchRatio(end)-1))) /...
    length(results.shuffledist_patchInterpatchRatio{end});