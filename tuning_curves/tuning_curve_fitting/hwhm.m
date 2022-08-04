%%%% HWHM: get right sided half width at half maximum above baseline iteratively
% arg1 = function handle of the form out=func(function_params,independent_variable)
% arg2 = function_params
% arg3 = min and max value of independent variable (x) params within which to search for hwhm
%        ...to gaurantee right-sided hwhm, min should be response peak val or higher 
% arg4 = spacing of independent variable points to check
% arg5 = independent variable value of peak response
% arg6 = baseline response rate
%%%% assumes downward curvature to right of response peak
%%%% if hwhm not found in range of x values specified, function returns NaN
%%%% updated 3/6/18 

function hwhm_out = hwhm(funchandle, funcparams, minmaxX, x_spacing, x_peak, baseline)


xvals = minmaxX(1):x_spacing:minmaxX(2);
yvals = funchandle(funcparams,xvals);
ypeak = funchandle(funcparams,x_peak);

ind = 2;
hwhmdone = false;
ysearch = yvals(1);
while (ysearch > (baseline + (ypeak-baseline)/2))
    ysearch = yvals(ind);
    if ind == length(yvals) % if hwhm hasn't been found before exhausting y values
        hwhm_out = NaN;
        hwhmdone = true;
        break;
    else
        ind = ind + 1;
    end
end

if ~hwhmdone  % if hwhm has been found before exhausting y values
    hwhm_out = xvals(ind-1) - x_peak;
end
    