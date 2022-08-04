%data from w030d (wally depth discrimination)
signed_x = [-64 -32 -16 -8 -4 4 8 16 32 64]'; %e.g., binoc corr

nostim_flags = [0 0 0 0 0 0 0 0 0 0]';
stim_flags = [1 1 1 1 1 1 1 1 1 1]';

nostim_slope = signed_x.*nostim_flags;
stim_slope = signed_x.*stim_flags;

pref_choice_nostim = [0 2 0 4 5 3 5 6 7 8]'; %# of preferred choices, no stimulation
pref_choice_stim = [0 3 5 8 6 9 7 8 10 8]'; %# of preferred choices, with stimulation

N_nostim = [10 10 10 10 10 10 10 10 10 10]';
N_stim = [10 10 11 10 10 10 10 10 10 10]';

X_nostim = [signed_x nostim_flags nostim_slope];
X_stim = [signed_x stim_flags stim_slope];
X = [X_nostim; X_stim];

Y_nostim = [pref_choice_nostim N_nostim];
Y_stim = [pref_choice_stim N_stim];
Y = [Y_nostim; Y_stim];

[b, dev, stats] = glmfit(X, Y, 'binomial', 'logit');

yfit_nostim = glmval(b, X_nostim, 'logit');
yfit_stim = glmval(b, X_stim, 'logit');

figure;

hold on;
plot(signed_x, pref_choice_nostim./N_nostim, 'bo');
plot(signed_x, yfit_nostim, 'b-');

plot(signed_x, pref_choice_stim./N_stim, 'ro');
plot(signed_x, yfit_stim, 'r-');
hold off;

line = sprintf('P corr = %10.8f', stats.p(2));
disp(line);
line = sprintf('P mstim = %10.8f', stats.p(3));
disp(line);
line = sprintf('P slope = %10.8f', stats.p(4));
disp(line);
