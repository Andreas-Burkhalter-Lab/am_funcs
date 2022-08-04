function[exptstats,trlstats,msacstats,intstats]=getmsacstats(msac,emt,emh,emv);

% function[exptstats,trlstats,msacstats,intstats]=getmsacstats(msac,emt,emh,emv);
% Function to take eye movement data and the output of the microsaccade-
%	picking program (microsaccade onset and offset times) and compute
%	statistics of the microsaccades for that experiment.
% INPUTS:  msac is zeros and ones; emt is times and nans; emh and emv are
%	eye position values in degrees.
% OUTPUTS:  outputs are organized into four matrices where each column
%	contains a parameter of interest and each row is an instance; each
%	output variable represents a different level of grouping --
%	whole experiment, by trial, by microsaccade, and by interval,
%	and each parameter of interest is placed in the appropriate variable.
%	-- exptstats contains stats for the whole experiment
%		col1	msaccount	number of saccades in experiment
%		col2	msacrate	rate of saccade occurrence
%	-- trlstats contains parameters for each trial
%		col1	trial number
%		col2	emtst		start of eye data for this trial
%		col3	emten		end of eye data for this trial
%	-- msacstats contains parameters for each microsaccade
%		col1	sacst		start time for saccade
%		col2	sacen		end time for saccade
%		col3	sacdur		duration
%		col4	sacamp		amplitude
%		col5	trial number
%		col6  sacdirec		direction of saccade
%		col7	sacspeed		average saccade speed
%	-- intstats contains parameters for each intersaccade interval
%		col1	imsi		intersaccade interval duration
%		col2	trial number
%		col3	rpflag		flag for real or pseudo interval
%					1 = real; 0 = pseudo
% Written by Crista, 10/18/98, modified by Greg and Ben 11/8/99

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Getting Variables %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Getting Trial Variables
tottrl=size(emt,1);
emtst=nanmin(emt')';
emten=nanmax(emt')';

%%% Getting Global Saccade Parameters
msaccount=length(find(msac==1));
tottime=sum(emten-emtst);
msacrate=msaccount/(tottime/1000);

%%% Getting Trial Numbers for Saccades
trlnums=[1:tottrl]'*ones(1,size(emt,2));
trlnums2=trlnums';
trlindex=trlnums2(find(msac'==1));
	
%%% Getting Individual Saccade Parameters
emt2=emt';
emh2=emh';
emv2=emv';
sacon=find(msac'==1);
sacoff=find(msac'==-1);

sacst=emt2(sacon);
sacen=emt2(sacoff);
sacdurs=(sacen-sacst);

sacsth=emh2(sacon);
sacstv=emv2(sacon);
sacenh=emh2(sacoff);
sacenv=emv2(sacoff);
sacamph=abs(sacenh-sacsth);
sacampv=abs(sacenv-sacstv);
sacamp=sqrt(sacamph.^2+sacampv.^2);
sacdirec = atan2((sacenv-sacstv), (sacenh-sacsth))*180.0/pi;

%average saccade speed in degrees/sec
avgsacspeed = 1000*(sacamp./sacdurs);

%%% Getting Intersaccade Interval Parameters
rimsis=[];
pimsis=[];
for trl=1:size(msac,1)
	msacst=emt(trl,find(msac(trl,:)==1));
	msacen=emt(trl,find(msac(trl,:)==-1));
	if length(msacst)>1
		subrimsis=msacst(2:length(msacst))-msacen(1:length(msacen)-1);
		subrimsis=[subrimsis' ones(length(subrimsis),1)*trl];
		rimsis=[rimsis;subrimsis];
	end
	if length(msacst)>0
		subpimsis=[msacst(1)-emtst(trl) emten(trl)-msacen(length(msacen))];
	else
		subpimsis=[emten(trl)-emtst(trl)];
	end
	subpimsis=[subpimsis' ones(length(subpimsis),1)*trl];
	pimsis=[pimsis;subpimsis];
end
rpflag=[ones(size(rimsis,1),1);zeros(size(pimsis,1),1)];
	
%%% Creating Output
exptstats=[msaccount msacrate];
trlstats=[[1:tottrl]' emtst emten];
msacstats=[sacst sacen sacdurs sacamp trlindex sacdirec avgsacspeed];
intstats=[[rimsis;pimsis] rpflag];
