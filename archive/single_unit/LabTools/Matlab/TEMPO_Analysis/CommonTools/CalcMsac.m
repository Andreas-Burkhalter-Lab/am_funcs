function[msac, emt]=CalcMsac(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol)

TEMPO_defs;
h = data.htb_header{EYE_DB};	%for convenience
eye_bin_width = 1000*(h.skip + 1) / (h.speed_units / h.speed);

emh= (squeeze(data.eye_data(REYE_H,:,BegTrial:EndTrial) ) )';
emv = (squeeze(data.eye_data(REYE_V,:,BegTrial:EndTrial) ) )';
emt = NaN* emh;

%% loop puts NaN outside of analysis window
for trial = 1:size(emh,1)
    time = StartEventBin(trial):eye_bin_width:StopEventBin(trial);
    index = floor(StartEventBin(trial)/eye_bin_width)  + 1;
    emt(trial,index: (index + length(time) - 1) ) = time;
end
% function[msac]=msacalg(emt,emh,emv,algparams,padding)
% Function applies an algorithm to eye movement data to determine
%	saccade onset and offset times.
% INPUTS:  emh, emv and emt are eyes position traces and times;
%	these variables are matrices of the same size:  row = trial
%	and columns are filled with consecutive data points;
%	algparams specifies the parameters of the algorithm:
%		algparams{1} = thrhigh     scalar, mandatory; deg/sec
%		algparams{2} = thrlow      scalar; deg/sec
%		algparams{3} = mindur      scalar; msec
%		algparams{4} = minisi      scalar; msec
%		algparams{5} = convfilt1   vector
%		algparams{6} = convfilt2   vector, mandatory
%	padding is an optional input that tells program what padding
%	value to expect at the end of the rows (e.g. nan, 0), however,
%	program will malfunction if expected padding value is incorrect.
%	To NOT specify an input or a cell in algparams, leave empty.
%	Defaults:  padding -- 0.
% ALGORITHM:  
%  	Horizontal (emh) and vertical (emv) eye positions are 
%	convolved with convfilt1 for smoothing purposes (to get 
%	femh and femv).  These, along with the time vector (emt), 
%	are convolved with convfilt2 to achieve difference values
%	(demh, demv, demt); convfilt2 should be a differentiator.
%	Vectorial velocity (vemt) is derived as follows:
%		vemt = sqrt(demh.^2 + demv.^2)/demt
%	Saccades are detected with thrhigh (high threshold in
%	degrees per second), and precise onset and offset times by
%	backing up or moving forward to thrlow (low threshold).  This
%	result is modified by eliminated detected saccades that are
%	under mindur (minimum duration in msec); consecutive saccades
%	with less than minisi (minimum intersaccade interval in msec)
%	are joined into one.
% OUTPUTS:  msac is a matrix of the same size as the emt with
%	a 1 at each saccade onset, -1 at each offset (actually
%	at the first timepoint after offset), and zeros elsewhere.
% Written by Crista, 8/18/98. Revised for TEMPO by BJP 7/20/01

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Checking Inputs %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaultalgparams=[{8};{6};{12};{20};{[.25 .5 .25]};{[-1 0 1]}];
defaultpadding=NaN;

algparams = defaultalgparams;
padding = defaultpadding;
%if nargin<4;algparams=defaultalgparams;end
%if isempty(algparams);algparams=defaultalgparams;end
%if nargin<5;padding=defaultpadding;end
%if isempty(padding);padding=defaultpadding;end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Setting Variables %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

femh=[];femv=[];
demh=[];demv=[];demt=[];
vemt=[];absvemr=[];
msac=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Getting Variables %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
trltot=size(emt,1);
samptot=size(emt,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Error Checking %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if (isempty(algparams{1}) | isempty(algparams{6}))
	disp('error in msacalg:  algparams field 1 or 6 is empty')
	algparams=[];
end

if ~isempty(algparams)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Getting Algorithm Paramters %%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	thrhigh=algparams{1}/1000;	% converting from deg/sec to deg/msec
	thrlow=algparams{2}/1000;	% converting from deg/sec to deg/msec
	mindur=algparams{3};		% in msec
	minisi=algparams{4};		% in msec
	convfilt1=algparams{5};
	convfilt2=algparams{6};
	
	if isempty(mindur);mindur=0;end
	if isempty(minisi);minisi=0;end
	
	if thrlow>thrhigh
		disp('warning in msacalg:  low threshold > high threshold')
		disp('   low threshold will be set to high threshold')
		thrlow=thrhigh;
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Converting Padding to Nans %%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if ~isnan(padding)
		x=find(emt==padding);
		emt(x)=nan;
	end
	
	x=find(isnan(emt));
	emh(x)=nan;
	emv(x)=nan;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Filtering Eye Position %%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if ~isempty(convfilt1)
			
		femh=conv2([1],[convfilt1],emh,'valid');
		femv=conv2([1],[convfilt1],emv,'valid');
		offfront=ceil((length(convfilt1)-1)/2);
		offback=floor((length(convfilt1)-1)/2);
		femh=[ones(size(femh,1),offfront)*nan femh ones(size(femh,1),offback)*nan];
		femv=[ones(size(femv,1),offfront)*nan femv ones(size(femv,1),offback)*nan];
		
	else

		femh=emh;
		femv=emv;

	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Calculating Difference Values %%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
	demh=conv2([1],[convfilt2],femh,'valid');
	demv=conv2([1],[convfilt2],femv,'valid');
	demt=conv2([1],[convfilt2],emt,'valid');
	offfront=ceil((length(convfilt2)-1)/2);
   offback=floor((length(convfilt2)-1)/2);
   
   %amplitudes in x and y direction, duration of movement - BJP
	demh=[ones(size(demh,1),offfront)*nan demh ones(size(demh,1),offback)*nan];
	demv=[ones(size(demv,1),offfront)*nan demv ones(size(demv,1),offback)*nan];
	demt=[ones(size(demt,1),offfront)*nan demt ones(size(demt,1),offback)*nan];
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Calculating Velocity %%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
   %calculate instantaneous velocities by taking square root of component amplitudes and dividing by the time between samples - BJP
   vemr=(sqrt(demh.^2+demv.^2))./demt;
   
   %convert to speeds - BJP
	absvemr=abs(vemr);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Applying High Threshold %%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %next line places 1 where eye speed is greater than thrhigh - BJP
	msaconh=absvemr>=thrhigh;
   
   % next line places a 1 or -1 in msac by calculating difference of indices in msaconh (1 = start of saccade, -1 = end of saccade) - BJP
	msac=[zeros(trltot,1) diff(msaconh')'];
	x=find(isnan(absvemr));msac(x)=0;
	for trl=1:trltot
      % check real indices after filtering of eye positions - BJP
      realindex=find(~isnan(absvemr(trl,:)));
	   realindexend=realindex(length(realindex));
      %don't want saccade to start at last index - BJP
      if msac(trl,realindexend)==1
			msac(trl,realindexend)=0;
      end
      %truncate saccade if it continues to end of index - BJP
		if msaconh(trl,realindexend)==1 & msaconh(trl,realindexend-1)==1
			msac(trl,realindexend)=-1;
		end
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Applying Low Threshold %%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
   if ~isempty(thrlow)
		msaconl=absvemr>=thrlow;
		for trl=1:trltot
		
			%%% Backing Up Starting Points
			msacst=find(msac(trl,:)==1);

			if ~isempty(msacst)			% if there were saccades
            
            % Moving Backward Starting points
            % Goal of next loop is to extend saccade before original starting point (obtained using thrhigh) 
            % to points with speed < thrhigh if thrlow is less than thrhigh - BJP
				msacst=msacst(msacst~=1);
				for col=msacst
					newcols=find(msaconl(trl,1:col)==0);
					if isempty(newcols)
						newcol=1;
					else
						newcol=newcols(length(newcols))+1;
					end
					msac(trl,newcol:col)=[1 zeros(1,col-newcol)];
				end
            
            % Goal of next loop is to extend saccade beyond original end point (obtained using thrhigh) 
            % to points with speed < thrlow if thrlow is less than thrhigh - BJP
            
				%%% Moving Forward End Points
				msacen=find(msac(trl,:)==-1);
				msacen=msacen(msacen~=samptot);
				for col=msacen
					newcols=find(msaconl(trl,col:samptot)==0)+col-1;
					if isempty(newcols)
						newcol=samptot;
					else
						if newcols(1)==col
							newcol=col;
						else
							newcol=newcols(1)-1;
						end
					end
					msac(trl,col:newcol)=[zeros(1,newcol-col) -1];
				end
				
			end
		end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Applying Minimum InterSaccade Interval %%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if minisi>0
		for trl=1:trltot
			msacen=find(msac(trl,:)==-1);
			if ~isempty(msacen)
				msacen=msacen(msacen~=samptot);
				if ~isempty(msacen)
					for col=msacen
						reltms=emt(trl,col+1:samptot)-emt(trl,col);
						reltms=find(reltms<minisi);
						if ~isempty(reltms)
							reltms=reltms(length(reltms));
							msacst=find(msac(trl,col:col+reltms)==1);
							if ~isempty(msacst)
								msacst=msacst(length(msacst))-1;
								msac(trl,col:col+msacst)=0;
							end
						end
					end
				end
			end
		end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Applying Minimum Duration %%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if mindur>0
		for trl=1:trltot
			msacst=find(msac(trl,:)==1);
			msacen=find(msac(trl,:)==-1);
			if ~isempty(msacst)
				sacdurs=emt(trl,msacen)-emt(trl,msacst);
				index=find((sacdurs<mindur) | (sacdurs>60)); %filter out blink artifacts with durations longer than 60 ms
				if ~isempty(index)
					msac(trl,msacst(index))=0;
					msac(trl,msacen(index))=0;
            end
			end
		end
   end   
end