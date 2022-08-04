% CHI2_TEST.M: Implements a chi-square goodness of fit test.  
%   datax, datay are the data to be fit (one x,y pair per trial)
%   funcname is the function that was fit, and params are the fit parameters to be tested
%   num_free_params is added so that you can correctly assess chi-square for functions in which there are
%   N free parameters, but some are held constant
%   GCD, starting 10/26/01
%function [chi2, chiP] = Chi2_Test(datax, datay, funcname, params, num_free_params)
function [chi2, chiP] = Chi2_Test_Fval(datax, datay, fity, num_free_params) %data already fitted

%make sure that datax, datay are column vectors
if ( (size(datax,1) == 1) & (size(datax,2) > 1) )   % if row vector
    datax = datax';
    datay = datay';
end

%find unique values of datax
uniquex = munique(datax);

%compute function values at uniquex
%tempaz = azimuth
%funccall = [funcname '(params, tempaz, tempel)'];
%fity = eval(funccall);
fity(fity < 0) = 0;

%to homogenize variance, work with the square roots of the data and fit
datay = sqrt(datay);
fity = sqrt(fity);

%compute the chi-square metric based on mean responses (weighted by s.d.)
std_err = []; std_dev = []; temp = []; 
sd_zero_count = 0;
for i=1:length(uniquex)
    select = logical(datax == uniquex(i));   
    std_err(i) = std(datay(select))/sqrt(length(select));
    std_dev(i) = std(datay(select));
    if (std_dev(i) == 0)
        sd_zero_count = sd_zero_count + 1;
        temp(i) = NaN;
    else
        temp(i) = (mean(datay(select)) - fity(i))^2 / std_dev(i)^2;
    end
end
chi2 = nansum(temp);
if (sd_zero_count > 0)
    str = sprintf('WARNING: %d values have been omitted from chi-sq computation due to zero standard deviation', sd_zero_count);
    disp(str);
end

%now, test significance of chi2
df_mean = (length(uniquex) - sd_zero_count) - num_free_params;  %degrees of freedom
    %Note, we subtract off sd_zero_count here to properly account for any omitted values of uniquex
    %due to zero std dev

if (df_mean <= 1)
    chiP = NaN;
    str2 = sprintf('CHI2_TEST ERROR: degrees of freedom (%d) not at least 1', df_mean);
    disp(str2);
else
    chiP = 1 - chi2cdf(chi2, df_mean);
end

return;