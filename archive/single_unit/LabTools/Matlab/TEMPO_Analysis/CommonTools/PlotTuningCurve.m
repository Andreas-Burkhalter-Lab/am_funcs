%----------------------------------------------------------------------------------------------------------
%-- PlotTuningCurve.m: This function takes a list of (x,y) pairs, computes the mean and s.e. of y at each
%-- 	unique value of x, and then plots the mean+/-s.e. as a plot with error bars.  
%--	Note that there are no calls to figure() here, so that this function can
%-- 	be invoked multiple times from a calling function as needed.	GCD, 1/3/2000
%--	NOTES:
%--		- function returns the x, y, and error values in vectors
%-- 		- calling the function with plot_flag=0 can be used to just get the x,y, and err values without graphing
%--		- I've added two output structs spk_max and spk_min, which correspond to the disp values and spike rates 
%--		  for max and min for either the spline fit or the data points if spline == [] - BJP 2/1/00
%----------------------------------------------------------------------------------------------------------
function	[unique_x, mean_rate, std_err, spk_max, spk_min] = PlotTuningCurve(x_list, y_list, symbol, line, spline, plot_flag);

mean_rate = []; std_err = [];

unique_x = munique(x_list);	%get unique values of x
for i=1:length(unique_x)
    matches = (x_list == unique_x(i));
    mean_rate(i) = mean( y_list(matches) );
    std_err(i) = std( y_list(matches) ) / sqrt(sum(matches));	%s.e. = sd/sqrt(N)
end
mean_rate = mean_rate';
std_err = std_err';

%plot the symbols and errorbars
if (plot_flag)
    errorbar(unique_x, mean_rate, std_err, std_err, symbol);
end

%if desired, plot lines and or splines
if (~isempty(line))
    if (spline)
        %interpolate the data using a cubic spline, then plot   
        x_interp = unique_x(1):0.01*(unique_x(length(unique_x)) - unique_x(1)): unique_x(length(unique_x));
        y_interp = interp1(unique_x, mean_rate, x_interp, 'spline');  %spline fit
        

        if (plot_flag)
            hold on;
            plot(x_interp, y_interp, line);
            hold off;
        end

    else
        if (plot_flag)
            hold on;
            plot(unique_x, mean_rate, line);
            hold off;
        end
    end
end

%find max and min of spline fit or of data points - BJP 2/1/00
if ((spline) & (~(isempty(line))))	%if fit curve to spline data
    [spk_max.y max_i] = max(y_interp);
    [spk_min.y min_i] = min(y_interp);
    spk_max.x = x_interp(max_i);
    spk_min.x = x_interp(min_i);
else    		% if fit no spline
    [spk_max.y max_i] = max(mean_rate);
    [spk_min.y min_i] = min(mean_rate);
    spk_max.x = unique_x(max_i);
    spk_min.x = unique_x(min_i);
end   


return;