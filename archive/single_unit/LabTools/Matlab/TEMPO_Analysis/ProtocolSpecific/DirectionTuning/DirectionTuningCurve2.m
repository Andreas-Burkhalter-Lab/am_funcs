%-----------------------------------------------------------------------------------------------------------------------
%-- DirectionTuningCurve2.m -- Plots a direction tuning curve and computes/plots a Gaussian fit to this curve
%--	GCD, 1/23/00
%-----------------------------------------------------------------------------------------------------------------------
function DirectionTuningCurve2(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

% not implemented yet to select output
output = 1;

ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

symbols = {'bo' 'r*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'b-' 'r-.' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};

% direction = data.obj_params(OBJ_ROTATION,:,1) .* data.obj_params(OBJ_DIR,:,1);
direction = data.obj_params(OBJ_ROTATION,:,1);


%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

    %some stuff done to send data to Izumi: GCD, 2/21/01
    %disp('saving...');
    %spikes_mat = squeeze(data.spike_data(1,:,:))';
    %size(spikes_mat)
    %outfid = fopen('temp2.dat', 'w');
    %for i=1:size(spikes_mat, 1)
    %    fprintf(outfid, '%d ', find(spikes_mat(i,:)>0) );
    %    fprintf(outfid, '\n');
    %end
    %fclose(outfid);


%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (data.misc_params(CONDITION,:) == 0) );  		%null condition = condition 0
%null_trials = logical( (direction == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(direction);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%   [direction' null_trials' select_trials']

dir(1,:) = data.obj_params(OBJ_DIR,:,1) == 1;
dir(2,:) = data.obj_params(OBJ_DIR,:,1) == -1;


figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 100 500 573], 'Name', 'Direction Tuning Curve');
subplot(2, 1, 2);

for set = 1:2

	plot_x = direction(~null_trials & select_trials & dir(set, :) );
	plot_y = spike_rates(~null_trials & select_trials &  dir(set, :) );


	%NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
	[px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{set}, lines{set}, 1, 0);	%note: last arg=0 means just get output, no plot

	unique_dirs = px; % store category groups for ANOVAs

	%now, shift the px, py, and perr vectors such that the peak of tuning curve is in middle of axis range
	% now, we need to shift these vectors so that the peak response is always in middle of vector
	ctr_indx = round(length(px)/2 - rem(length(px),2)) + 1;
	[max_val max_indx] = max(py);
   shift = max_indx - ctr_indx;
 
   if (shift > 0)
   	px = [ px(shift+1 : length(px)) ; px(1 : shift)+180];
    	py = [ py(shift+1 : length(py)) ; py(1 : shift)];
      perr = [ perr(shift+1 : length(perr)) ; perr(1 : shift)];
	end
	if (shift < 0)
   	px = [ px(length(px)+shift+1 : length(px))-180 ; px(1 : length(px)+shift)];
    	py = [ py(length(py)+shift+1 : length(py)) ;  py(1 : length(py)+shift)];
    	perr = [ perr(length(perr)+shift+1 : length(perr)) ;  perr(1 : length(perr)+shift)];
   end


	% Since direction is circular, duplicate  the lowest value (px) as px+ 360, and change spike arrays accordingly
	% px = [px; px(1)+360];
	% py = [py; py(1)];
	% perr = [perr; perr(1)];

	hold on;
	%plot the data, after shifting as necessary above
	errorbar(px, py, perr, perr, symbols{set} );
	hold on;

	%now, fit the data with a Gaussian curve and plot this as well
	means = [px py];
	raw = [px py];

	[pars] = gaussfit(means, raw, 1);
	x_interp = (px(1)): 1 : (px(length(px)));
	y_interp = gaussfunc(x_interp, pars);
    yl(set) = max(y_interp');
   
	h(set) = plot(x_interp, y_interp, lines{set});

	%now, get the firing rate for NULL condition trials and add spontaneous rate to plot
	null_x = [min(px) max(px)];
	null_rate = mean(data.spike_rates(SpikeChan, (null_trials & select_trials)  ));
	null_y = [null_rate null_rate];
	hold on;
	plot(null_x, null_y, 'k--');
	hold off;



	% calculate some metrics and stats then print them in plot
	pref_dir(set) = pars(3);
	p_value(set) = spk_anova(plot_y, plot_x, unique_dirs);
	base_rate(set) = pars(1);
	amplitude(set) = pars(2);
	max_rate(set) = base_rate(set) + amplitude(set);
	width(set) = sqrt(-(log(.5)))*pars(4)*2*sqrt(2);
	DSI(set) = 1 - (base_rate(set) - null_rate)/(max_rate(set) - null_rate); 

end

yl = max(yl);
YLim([0 yl]);	% set the lower limit of the Y axis to zero
XLabel('Direction of Motion (deg)');
YLabel('Response (spikes/sec)');


%now, print out some useful information in the upper subplot
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

if max_rate(1) > max_rate(2)
   PrintDirectionData(p_value(1), base_rate(1), null_rate, amplitude(1), pref_dir(1), max_rate(1), width(1), DSI(1)); 
else 
   PrintDirectionData(p_value(2), base_rate(2), null_rate, amplitude(2), (pref_dir(2) + 180) , max_rate(2), width(2), DSI(2)); 
end


legend(h,'0<dir<180','180<dir<360', 4);

return;