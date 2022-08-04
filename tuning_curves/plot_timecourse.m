




replace timecourse mean with 'pref_timecourse' variable, which needs to be stored in tuning table


timecourse_mean = [tuningdat.tuning.sf_trials{4}.timecourse_mean{3,1}, tuningdat.tuning.sf_trials{4}.timecourse_mean{3,2}, tuningdat.tuning.sf_trials{4}.timecourse_mean{3,3}]

plot(timecourse_mean)
hold on
errorbar(1:length(timecourse_mean),timecourse_mean,...
    [tuningdat.tuning.sf_trials{4}.timecourse_sem{3,1}, tuningdat.tuning.sf_trials{4}.timecourse_sem{3,2}, tuningdat.tuning.sf_trials{4}.timecourse_sem{3,3}])