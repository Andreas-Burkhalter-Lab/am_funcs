function[timsampsindex,ERROR]=regtimsamp(timsamps,testper,mode,padding,dispmode)

% function[timsampsindex,ERROR]=regtimsamp(timsamps,testper,mode,padding,dispmode)
% Function to regularize a set of time samples.  Searches for the 
%	first offset for which there are regular samples at the
%	sampling period for the whole trial (mode 1) or the longest
%	stretch of samples that exist at the sampling period of
%	all the offsets (mode 2).
% INPUTS:  timsamps is a matrix of time samples (in msec) where each
%	row = a trial and each column contains consecutive time
%	samples, padded with zeros at the end.  Testper is an optional
%	sampling period input; default = 4 if empty or not present.
%	Mode is an optional mode input (see above); default = 1.
%	Padding is value time matrix is padded with; default = 0.
%	Dispmode controls trial by trial display; 0=no, 1=yes.
% OUTPUTS:  timsampsindex is a matrix of the same size as timsamps
%	with ones for valid samples and zeros (or padding value if it
%	is not zero) elsewhere.  If any error messages are received, 
%	these samples may not be regular.  ERROR=1 means lack of regular
%	samples; ERROR=0 means success.
% Written by Crista, 8/28/98.

ERROR=0;

%%%%% Checking Inputs %%%%%
defaultmode=1;						% default mode
defaulttestper=4;					% default test period = 4 msec
defaultpadding=0;					% default padding of time matrix = 0
defaultdispmode=1;					% default dispmode is to display
if nargin<2;testper=defaulttestper;end
if isempty(testper);testper=defaulttestper;end
if nargin<3;mode=defaultmode;end
if isempty(mode);mode=defaultmode;end
if nargin<4;padding=defaultpadding;end
if isempty(padding);padding=defaultpadding;end
if nargin<5;dispmode=defaultdispmode;end
if isempty(dispmode);dispmode=defaultdispmode;end

%%%%% Error Checking %%%%%
if isempty(timsamps)
	ERROR=1;
end

%%%%% Getting and Setting Variables %%%%%
trltot=size(timsamps,1);				% number of trials
maxnumsamp=size(timsamps,2);				% maximum number of samples
timsampsindex=zeros(trltot,maxnumsamp);			% output matrix

%%%%% Looping Through Trials %%%%%
for trl=1:trltot
	if isnan(padding)
		realtimsamps=find(~isnan(timsamps(trl,:)));
		timsamp=timsamps(trl,realtimsamps);
	else
		realtimsamps=find(timsamps(trl,:)~=padding);
		timsamp=timsamps(trl,realtimsamps);	% get time sample values for this trial
	end
	minsamp=min(timsamp);				% minimum time value for this trial
	maxsamp=max(timsamp);				% maximum time value for this trial
	lensamp=length(timsamp);			% number of time values in this trial

	%%%%% Searching for Regular Samples for Whole Trial %%%%%
	if mode==1
	
		%%% set offset found flag to 0
		perfound=0;
		
		%%% for each possible offset from the first time sample (minsamp)
		for sampstart=1:testper
		
			%%% create vector of expected time samples if there is a match
			testsamp=[minsamp+sampstart-1:testper:maxsamp];
			
			%%% vector of flags for which time values match
			check=ismember(testsamp,timsamp);
			
			%%% if all times match
			if length(find(check))==length(testsamp)
			
				%%% if this is the first offset found, store it in output
				if perfound==0
					backcheck=ismember(timsamp,testsamp);
					timsampsindex(trl,1:length(backcheck))=backcheck;
					if dispmode
						disp(['   trial ',num2str(trl),':  sample period ', ....
						num2str(testper),' at offset ',num2str(sampstart)]);
					end
								
				%%% if this is not the first offset found, display warning
				elseif perfound==1
					if dispmode
						disp(['   warning:  sample period ',num2str(testper), ....
						' also found at offset ',num2str(sampstart)])
					end
				end
				
				%%% set offset found flag to 1
				perfound=1;
			end
		end

		%%% if no offset was found for this trial
		if perfound==0
			if dispmode
				disp(['   trial ',num2str(trl),':  failed'])
			end
			ERROR=1;
		end
	
	%%%%% Searching for Longest Stretch of Regular Samples %%%%%
	elseif mode==2
	
		%%% sample search results matrix
		%%% columns represent each possible offset
		%%% rows contain maximum time left at this offset (timeleft),
		%%%   time lost (timelost), and sample start and end times
		%%%   (testsampstart and testsampend)
		%%% once this matrix is full, the maximum time out of all
		%%%   possible offset will be found and used
		ssres=zeros(4,testper);
		
		%%% for each possible offset from the first time sample (minsamp)
		for sampstart=1:testper
		
			%%% create vector of expected time samples if there is a match
			testsamp=[minsamp+sampstart-1:testper:maxsamp];
			
			%%% vector of flags for which time values match
			check=ismember(testsamp,timsamp);
			
			%%% vector of flags for start and end of continuous samples
			dcheck=[0 diff(check)];
			x=find(dcheck==1);			% start times
			y=find(dcheck==-1);			% end times
			
			%%% if there is a start and/or end time within time sample
			if ~isempty(x) | ~isempty(y)
			
				%%% checks that deal with ends of time sample
				if isempty(y)
					y=length(dcheck);
				end
				if isempty(x)
					x=[1];
				end
				if y(length(y))<x(length(x))
					y=[y length(dcheck)];
				end
				if y(1)<x(1);x=[1 x];end
				
				%%% vector of length in msec of continous samples within time sample
				z=y-x;
				
				%%% find the largest continous sample
				index=find(z==max(z));index=index(1);
				
				%%% get start and end times for this sample, put in ssres
				x=x(index);y=y(index);
				testsampstart=testsamp(x);
				testsampend=testsamp(y);
				timeleft=testsampend-testsampstart;
				timelost=(maxsamp-testsampend)+(testsampstart-minsamp);
				ssres(:,sampstart)=[timeleft;timelost;testsampstart;testsampend];

			%%% if there are no start or end times within time sample
			else
			
				%%% check if all samples are present; if so take whole sample
				if isempty(find(check~=1))
					testsampstart=testsamp(1);
					testsampend=testsamp(length(testsamp));
					timeleft=testsampend-testsampstart;
					timelost=(maxsamp-testsampend)+(testsampstart-minsamp);
					ssres(:,sampstart)=[timeleft;timelost;testsampstart;testsampend];
				end		
			end
		end
		
		%%% find largest continuous sample out of all possible offsets
		sampfind=find(ssres(1,:)==max(ssres(1,:)));sampfind=sampfind(1);
		
		%%% take that time sample
		testsampstart=ssres(3,sampfind);
		testsampend=ssres(4,sampfind);
		timeleft=ssres(1,sampfind);
		timelost=ssres(2,sampfind);

		%%% write this sample into output matrix
		testsamp=[testsampstart:testper:testsampend];
		backcheck=ismember(timsamp,testsamp);
		timsampsindex(trl,realtimsamps)=backcheck;
		
		if dispmode
			disp(['trl ',num2str(trl),'  min/maxsamp ',num2str(minsamp),' ',num2str(maxsamp), ....
			'  sampst/end ',num2str(testsampstart),' ',num2str(testsampend),'  timeleft/lost ', ....
			num2str(timeleft),' ',num2str(timelost)])
		end
	end

end

%%%%% Filling In Output With Padding Value %%%%%
if padding~=0
	dum=find(timsampsindex==0);
	timsampsindex(dum)=padding;
end

%%%%% Checking Results %%%%%
newtimsamps=timsamps.*timsampsindex;
for trl=1:trltot
	if isnan(padding)
		timsamp=newtimsamps(trl,find(~isnan(newtimsamps(trl,:))));
	else
		timsamp=newtimsamps(trl,find(newtimsamps(trl,:)~=padding));	% get time sample values for this trial
	end
	if isempty(timsamp)
		ERROR=1;
		if dispmode
			disp(['regtimsamp error:  trl ',num2str(trl),' is empty'])
		end
	else
		if std(diff(timsamp))~=0
			ERROR=1;
			if dispmode
				disp(['regtimsamp error:  trl ',num2str(trl),' is not regular'])
			end
		end
	end
end
