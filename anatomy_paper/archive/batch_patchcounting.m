%%%% batch patchcounting
% use plot_patchcount instead if images do not cover full extent of area
%   being processed 

function res = batch_patchcounting(intable,arealabel)
% rf mean diameter (deg) from RFs by areas from RD and WJ, jan 2018	
% por	43
% v1	8
% lm	12.5
% li	32

rf_diam_deg_por = 43; % por 
rf_diam_deg_v1 = 8; % v1
rf_diam_deg_lm = 12.5; % lm
rf_diam_deg_li = 32; % li

 % use a single roi area value to calculate patches per rf, rather than individual roi areas for each case
use_mean_roi_area = 0;

visualfield_sqdeg = 90*105; % approximate from wang and burkhalter 2007 fig10b


switch arealabel
    case 'por'
        rf_diam_deg = rf_diam_deg_por;
    case 'v1'
        rf_diam_deg = rf_diam_deg_v1;
    case 'lm'   
        rf_diam_deg = rf_diam_deg_lm;
    case 'li'
        rf_diam_deg = rf_diam_deg_li;
end

nfiles = height(intable);

for i = 1:nfiles
%     intable.patches{i} = analyzePatches(intable.m2{i},intable.roi{i},intable.zoom(i),intable.scope{i});
    intable.roi_area_sqmm(i,1) = intable.patches{i}.roi_area_sqmm;
    intable.patches_per_sqmm(i,1) = intable.patches{i}.patches_per_sqmm;
    if use_mean_roi_area % average roi size from all cases
        sqmm_per_sqdeg = mean(intable.roi_area_sqmm) / visualfield_sqdeg;
    else % case-specific roi size
        sqmm_per_sqdeg = intable.roi_area_sqmm(i) / visualfield_sqdeg; 
   
    end
    rf_area_sqdeg = pi * (rf_diam_deg/2)^2;
    rf_area_sqmm = rf_area_sqdeg * sqmm_per_sqdeg; % point image area
    rfs_per_sqmm = 1/rf_area_sqmm;
    intable.patches_per_rf(i) = intable.patches_per_sqmm(i) / rfs_per_sqmm;
    
end

res = intable;