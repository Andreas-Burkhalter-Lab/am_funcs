% checkgregdsdata.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% changing into the directory where the functions are %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fxnpath='/HOME/crista/matlab/eyes';
eval(['cd ',fxnpath])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% defining filename and paths for radar's data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename='rd043a';
datapath='/d5/gregd/depth/radar/spikes';
indpath='/HOME/gregd/Matlab';
algsetpath='/HOME/crista/matlab/eyes/algparamsets';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% defining other default parameters for radar's data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

regflg=2;		% flag to regularize the data
			% 0 = don't regularize; 1 = regularize the whole
			%	trial; 2 = regularize longest stretch of
			%	trial that you can find
truncatecol=4;		% column of index file to use in truncating
			% data at a particular event
			% column 4 = 'SIN_STOFF' for paradigm 364
cullcol=4;		% column of foo file to use to select subset of trials
			% column 4 = 'I_CTR_SPEED' for paradigm 364
cullmin=0;		% minimum value to use in selecting subset of trials
			% based on values in <cullcol> column of foo file
cullmax=80;		% maximum value to use in selecting subset of trials
			% based on values in <cullcol> column of foo file
			
algparams=cell(6,1);
algparams{1}=8;			% upper (initial) threshold velocity in deg/sec
algparams{2}=6;			% lower (second) threshold velocity in deg/sec
algparams{3}=12;		% minimum duration of microsaccade
algparams{4}=20;		% minimum interval between microsaccades (otherwise
				%	saccades are joined into one)
algparams{5}=[.25 .5 .25];	% filter on horz and vert eye position
algparams{6}=[-1 0 1];		% convolution filter used on emt, emh and emv

%%%%%%%%%%%%%%%%%%%%
%%% getting data %%%
%%%%%%%%%%%%%%%%%%%%

[foo,datindex,raster]=load_FooIndexRaster(filename,datapath);
eval(['load ',datapath,'/',filename,'.tim'])
eval(['emt=',filename,';'])
eval(['load ',datapath,'/',filename,'.e1h'])
eval(['emh=',filename,'./40;'])
eval(['load ',datapath,'/',filename,'.e1v'])
eval(['emv=',filename,'./40;'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% displaying raw data sizes %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp(['Raw Data Sizes:'])
disp(['        rows   cols'])
disp(['foo     ',num2str(size(foo,1)),'     ',num2str(size(foo,2))])
disp(['index   ',num2str(size(datindex,1)),'     ',num2str(size(datindex,2))])
disp(['raster  ',num2str(size(raster,1)),'     ',num2str(size(raster,2))])
disp(['emt     ',num2str(size(emt,1)),'     ',num2str(size(emt,2))])
disp(['emh     ',num2str(size(emh,1)),'     ',num2str(size(emh,2))])
disp(['emv     ',num2str(size(emv,1)),'     ',num2str(size(emv,2))])
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%
%%% processing data %%%
%%%%%%%%%%%%%%%%%%%%%%%

%%% (1)  removing redundant data points from eye movement data %%%

		%%%%% Removing Redundant Data Points %%%%%
		trltot=size(emt,1);
		newemt=ones(trltot,size(emt,2))*nan;
		newemh=ones(trltot,size(emt,2))*nan;
		newemv=ones(trltot,size(emt,2))*nan;
		for trl=1:trltot
			emtdat=emt(trl,emt(trl,:)~=0);
			emhdat=emh(trl,emt(trl,:)~=0);
			emvdat=emv(trl,emt(trl,:)~=0);
			t=[1 diff(emtdat)];
			tkeep=find(t~=0);
			emtdat=emtdat(tkeep);
			emhdat=emhdat(tkeep);
			emvdat=emvdat(tkeep);
			newemt(trl,1:length(emtdat))=emtdat;
			newemh(trl,1:length(emtdat))=emhdat;
			newemv(trl,1:length(emtdat))=emvdat;
		end
		clear emhdat emvdat emtdat t tkeep trltot trl
		
		%%%%% Removing Extra Zeros %%%%%
		x=nansum(newemt);
		y=size(newemt,2);
		while isnan(x(y));y=y-1;end
		newemt=newemt(:,1:y);
		newemh=newemh(:,1:y);
		newemv=newemv(:,1:y);
		clear x y
	
		%%%%% Reassigning Variables %%%%%
		emt=newemt;clear newemt
		emh=newemh;clear newemh
		emv=newemv;clear newemv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% displaying processed data sizes %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp(['Processed Data Sizes (after removing redundant data points):'])
disp(['        rows   cols'])
disp(['emt     ',num2str(size(emt,1)),'     ',num2str(size(emt,2))])
disp(['emh     ',num2str(size(emh,1)),'     ',num2str(size(emh,2))])
disp(['emv     ',num2str(size(emv,1)),'     ',num2str(size(emv,2))])
disp(' ')

%%% (2) setting paradigm number %%%

paradigm=364;

%%% (3) getting fixation point x and y coordinates %%%

fpxstring='FP X-Coord:';
fpystring='FP Y-Coord:';
index=0;fpx=[];fpy=[];
fid=fopen([datapath,'/',filename,'.s']);   % opening .s file

		%%%%% Obtaining Fixation Point Position Values %%%%%
		index=0;
		while index==0
			line=fgets(fid);
			if line==-1;index=1;end
			if length(line)>length(fpxstring)
				index=strcmp(line(1:length(fpxstring)),fpxstring);
				if index==1
					fpx=str2num(line(length(fpxstring)+1:length(line)))/10;
					if paradigm==363;fpx=fpx*(-1);end		% kluge
					if paradigm==364;fpx=fpx*(-1);end		% kluge
				end
				index=0;
			end
			if length(line)>length(fpystring)
				index=strcmp(line(1:length(fpystring)),fpystring);
				if index==1
					fpy=str2num(line(length(fpystring)+1:length(line)))/10;
				end
				index=0;
			end
			if ~isempty(fpx) & ~isempty(fpy)
				index=1;
			end
		end

		%%%%% Normalizing Data %%%%%
		emh=emh-fpx;
		emv=emv-fpy;

		fclose(fid);clear fid
		clear index fpx fpy fpxstring fpystring line fpstrings

%%% (4) getting more heavily processed data %%%

pemt=emt;
pemh=emh;
pemv=emv;
msac=[];
trlnums=[];

%%% (4a) getting subset of trials

		%%%%% Selecting Subset of Trials %%%%%
		if cullcol~=0				% if subset of trials is desired
			trlnums=[1:size(emt,1)]';
			if isempty(cullmin);cullmin=min(foo(:,cullcol));end
			if isempty(cullmax);cullmax=max(foo(:,cullcol));end
			testfoo=find(foo(:,cullcol)>=cullmin & foo(:,cullcol)<=cullmax);
			trlnums=trlnums(testfoo);
			pemt=pemt(trlnums,:);
			pemh=pemh(trlnums,:);
			pemv=pemv(trlnums,:);
		else
			trlnums=[1:size(emt,1)]';
		end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% displaying processed data sizes %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp(['Processed Data Sizes (after selecting subset of trials):'])
disp(['        rows   cols'])
disp(['pemt    ',num2str(size(pemt,1)),'     ',num2str(size(pemt,2))])
disp(['pemh    ',num2str(size(pemh,1)),'     ',num2str(size(pemh,2))])
disp(['pemv    ',num2str(size(pemv,1)),'     ',num2str(size(pemv,2))])
disp(' ')

%%% (4b) truncating data at event times specified in <truncatecol> column of index file %%%

		%%%%% Truncating Data %%%%%
		if truncatecol~=0
	
			validvalues=datindex(trlnums,truncatecol);
			dumcheck=nanmax(nanmax(pemt));
			dum=find(validvalues<0);validvalues(dum)=dumcheck;
			validtimes=validvalues*ones(1,size(pemt,2));
			validtimes=pemt<=validtimes;
			dum=find(validtimes==0);validtimes(dum)=nan;
			pemt=pemt.*validtimes;			% multiply invalid times by nans
			pemh=pemh.*validtimes;
			pemv=pemv.*validtimes;
			clear validvalues dumcheck dum validtimes
			
		end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% displaying processed data sizes %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp(['Processed Data Sizes (after truncating data):'])
disp(['        rows   cols'])
disp(['pemt    ',num2str(size(pemt,1)),'     ',num2str(size(pemt,2))])
disp(['pemh    ',num2str(size(pemh,1)),'     ',num2str(size(pemh,2))])
disp(['pemv    ',num2str(size(pemv,1)),'     ',num2str(size(pemv,2))])
disp(' ')

%%% (4c) regularizing the eye movement data time samples (removing irregular samples) %%%

		%%%%% Regularizing Data %%%%%
		if regflg~=0

			testper=4;dispmode=0;padding=nan;
			[pemtindex,ERROR2]=regtimsamp(pemt,testper,regflg,padding,dispmode);
			if ERROR2==0
				pemt=pemt.*pemtindex;
				pemh=pemh.*pemtindex;
				pemv=pemv.*pemtindex;
			else
				disp('error in calculating:  unable to regularize samples')
			end
			clear testper dispmode pemtindex
	
			%%%%% Removing Redundant Data Points %%%%%
			trltot=size(pemt,1);
			newemt=ones(trltot,size(pemt,2))*nan;
			newemh=ones(trltot,size(pemt,2))*nan;
			newemv=ones(trltot,size(pemt,2))*nan;
			for trl=1:trltot
				emtdat=pemt(trl,~isnan(pemt(trl,:)));
				emhdat=pemh(trl,~isnan(pemt(trl,:)));
				emvdat=pemv(trl,~isnan(pemt(trl,:)));
				t=[1 diff(emtdat)];
				tkeep=find(t~=0);
				emtdat=emtdat(tkeep);
				emhdat=emhdat(tkeep);
				emvdat=emvdat(tkeep);
				newemt(trl,1:length(emtdat))=emtdat;
				newemh(trl,1:length(emtdat))=emhdat;
				newemv(trl,1:length(emtdat))=emvdat;
			end
			clear emhdat emvdat emtdat t tkeep trltot trl
			
			%%%%% Removing Extra Zeros %%%%%
			x=nansum(newemt);
			y=size(newemt,2);
			while isnan(x(y));y=y-1;end
			newemt=newemt(:,1:y);
			newemh=newemh(:,1:y);
			newemv=newemv(:,1:y);
			clear x y
		
			%%%%% Reassigning Variables %%%%%
			pemt=newemt;clear newemt
			pemh=newemh;clear newemh
			pemv=newemv;clear newemv
	
		end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% displaying processed data sizes %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp(['Processed Data Sizes (after regularizing data):'])
disp(['        rows   cols'])
disp(['pemt    ',num2str(size(pemt,1)),'     ',num2str(size(pemt,2))])
disp(['pemh    ',num2str(size(pemh,1)),'     ',num2str(size(pemh,2))])
disp(['pemv    ',num2str(size(pemv,1)),'     ',num2str(size(pemv,2))])
disp(' ')

%%% (4d) getting microsaccade times

		%%%%% Getting Microsaccades %%%%%
		padding=nan;
		msac=msacalg(pemt,pemh,pemv,algparams,padding);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% displaying processed data sizes %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp(['Processed Data Sizes (microsaccade matrix):'])
disp(['        rows   cols'])
disp(['msac    ',num2str(size(msac,1)),'     ',num2str(size(msac,2))])
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% getting saccade triggered average %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subraster=raster(trlnums,:);
[msacta,msactatimes,msaccount]=msactamaker(msac,pemt,subraster,filename);



