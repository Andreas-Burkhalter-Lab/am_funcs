%---------------------------------------------------------------------------------------------------------------------
%-- PrintGeneralData.m : THis function prints out some useful general information at the top of each figure plot.
%--	GCD, 1/26/00
%---------------------------------------------------------------------------------------------------------------------

% 3/14/16 AM fixed typos

function   PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		%needed for defns..... % 3/14/16 AM chaned from 'Tempo_Defs'
ProtocolDefs;	%references to protocol Names 1/4/01 BJP

axis([0 100 0 100]);
axis('off');
xpos = -10;
ypos = 110;
font_size = 8;
bump_size = 7;

temp = strcat(PATH, FILE);
temp(temp == '\') = '/';
% this prevents a stupid error from appearing on the screen
line = sprintf('File: %s', temp);
%Next line is the source of stupid error
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Collected on: %s', data.htb_header{EVENT_DB}.date);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('TEMPO program: %s', data.htb_header{EVENT_DB}.pro_file);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('TEMPO config. file: %s', data.htb_header{EVENT_DB}.cfg_file);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
if (~isempty(data.eye_data))
    eh = data.htb_header{EYE_DB};
    line = sprintf('Eye data: %d channel(s) @ %6.1f Hz', eh.nchannels, eh.speed_units/(eh.speed * (eh.skip + 1)) );
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
end
if (~isempty(data.spike_data))
    sh = data.htb_header{SPIKE_DB};
    line = sprintf('Spike data: %d channel(s) @ %6.1f Hz', sh.nchannels, sh.speed_units/(sh.speed * (sh.skip + 1)) );
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
end
evh = data.htb_header{EVENT_DB};
line = sprintf('Event data: %d channel(s) @ %6.1f Hz', evh.nchannels, evh.speed_units/(evh.speed * (evh.skip + 1)) );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Good Trials: %d out of %d total', size(data.event_data, 3), evh.sweep);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Trial range for analysis: %d -> %d', BegTrial, EndTrial);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Spike channel analyzed: %d', SpikeChan);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Analysis Start Event: %s, offset %d ms', event_names{StartCode}, StartOffset);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Analysis Stop Event: %s, offset %d ms', event_names{StopCode}, StopOffset);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

xpos = 55;
ypos = 110;
font_size = 8;

%protocol number and print out neuron params accordingly
if (size(data.spike_data,1) > 2) %multielectrode protocols
	line = sprintf('Protocol: %d (%s)', Protocol, protocol_names{Protocol+1});
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
	line = sprintf('Pref. Direction: %6.1f %6.1f (deg)', data.neuron_params(PREFERRED_DIRECTION, :));
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
	line = sprintf('Pref. Speed: %6.1f %6.1f (deg/sec)', data.neuron_params(PREFERRED_SPEED, :));
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
%	line = sprintf('Pref. Horiz. Disp.: %6.2f %6.2f (deg)', data.neuron_params(PREFERRED_HDISP, :));
%	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
	line = sprintf('Pref. Orientation: %6.1f %6.1f (deg)', data.neuron_params(PREFERRED_ORIENTATION, :));
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
	line = sprintf('Pref. SF: %6.3f %6.3f (cyc/deg)', data.neuron_params(PREFERRED_SPATIAL_FREQ, :));
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
	line = sprintf('Pref. TF: %6.2f %6.2f (Hz)', data.neuron_params(PREFERRED_TEMPORAL_FREQ, :));
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
	line = sprintf('RF X-Ctr (rel to FP): %6.1f %6.1f (deg)', data.neuron_params(RF_XCTR, :));
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
	line = sprintf('RF Y-Ctr (rel to FP): %6.1f %6.1f (deg)', data.neuron_params(RF_YCTR, :));
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
	line = sprintf('RF Diameter: %6.1f %6.1f (deg)', data.neuron_params(RF_DIAMETER, :));
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size
else
	line = sprintf('Protocol: %d (%s)', Protocol, protocol_names{Protocol+1});
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
	line = sprintf('Pref. Direction: %6.1f (deg)', data.neuron_params(PREFERRED_DIRECTION, 1));
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
	line = sprintf('Pref. Speed: %6.1f (deg/sec)', data.neuron_params(PREFERRED_SPEED, 1));
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
%	line = sprintf('Pref. Horiz. Disp.: %6.2f (deg)', data.neuron_params(PREFERRED_HDISP, 1));
%	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
	line = sprintf('Pref. Orientation: %6.1f (deg)', data.neuron_params(PREFERRED_ORIENTATION, 1));
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
	line = sprintf('Pref. SF: %6.3f (cyc/deg)', data.neuron_params(PREFERRED_SPATIAL_FREQ, 1));
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
	line = sprintf('Pref. TF: %6.2f (Hz)', data.neuron_params(PREFERRED_TEMPORAL_FREQ, 1));
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
	line = sprintf('RF X-Ctr (rel to FP): %6.1f (deg)', data.neuron_params(RF_XCTR, 1));
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
	line = sprintf('RF Y-Ctr (rel to FP): %6.1f (deg)', data.neuron_params(RF_YCTR, 1));
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
	line = sprintf('RF Diameter: %6.1f (deg)', data.neuron_params(RF_DIAMETER, 1));
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
end;

return;