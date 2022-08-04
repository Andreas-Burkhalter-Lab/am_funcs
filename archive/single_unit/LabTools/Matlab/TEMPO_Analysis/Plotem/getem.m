function[varout1,varout2,varout3,varout4,varout5,varout6,varout7]=getem(filename,outputlist)

% function[varout1,varout2,varout3,varout4,varout5,varout6,varout7]=getem(filename,outputlist)
% Function to get eye movement data.  Based on plotem.m but without
%	the graphics; allows user to specify output variables.
% INPUT:  filename is a string, e.g. 'ds014a'.; outputlist is 
%	a list of names from the list below specifying
%	which variables to output, and in what order.
% 	Output list:  foo, datindex, raster, emt, emh, emv, pemt,
%		pemh, pemv, msac, msacstats, msacta, msactatimes,
%		msaccount; indexcols, eventnames.
% LOADS:  <filename>.s, .index, .ras, .tim, .e1h, .e1v
% CALLS:  plotemdefaults.m, load_FooIndexRaster.m, regtimsamp.m,
%	  msacalg.m
% Written by Crista, 10/13/98.
% Datapath structure modified, 1/7/99.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Setting Default Variables %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ERROR=0;
varout1=[];varout2=[];varout3=[];varout4=[];varout5=[];varout6=[];varout7=[];
currpath=cd;
fxnpath='/HOME/crista/matlab/eyes';
eval(['cd ',fxnpath])

eventnames={'None';
	    'SIN_FPON';
	    'SIN_FIX';
	    'SIN_STON';
	    'SIN_STOFF';
	    'SIN_TRGON';
	    'SIN_TRGOFF';
	    'SIN_FPOFF';
	    'SIN_SAC';
	    'SIN_TRGAC';
	    'SIN_DTON';
	    'SIN_DTOFF';
	    'SIN_SACST';
	    'SIN_SACEN'};
varnames={'foo';
	  'index';
	  'raster';
	  'emt';
	  'emh';
	  'emv';
	  'msac';
	  'pemt';
	  'pemh';
	  'pemv';
	  'pmsac'};
fxncalls={['[summarystats]'],['plotmsacstats(pmsac,pemt,pemh,pemv,[],0)'];
	  ['[msacta,msactatimes,msaccount]'],['msactamaker(pmsac,pemt,raster,filename,0)']};


disp(['FILE:  ',filename])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Loading Default Parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('   loading defaults')

plotemdefaults

if (regflg<0 | regflg>2);regflg=0;end
truncatecol=1;
for x=1:size(eventnames,1)
	if strcmp(truncatestring,char(eventnames(x)))
		truncatecol=x;
	end
end
tzerovalue=0;
for x=1:size(eventnames,1)
	if strcmp(tzerostring,char(eventnames(x)))
		tzerovalue=x;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Loading Data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
disp('   loading data')

foo=[];datindex=[];raster=[];emt=[];emh=[];emv=[];
pemt=[];pemh=[];pemv=[];msac=[];
paradigm=[];fpxstring=[];fpystring=[];

for dpi=1:size(datapath,2)
  currdatapath = char(datapath(dpi));
  if (isempty(foo) & exist([currdatapath,'/',filename,'.s'])==2)
    [foo,datindex,raster]=load_FooIndexRaster(filename,currdatapath);
  end
  if (isempty(emt) & exist([currdatapath,'/',filename,'.tim'])==2)
    eval(['load ',currdatapath,'/',filename,'.tim'])
    eval(['emt=',filename,';'])
  end
  if (isempty(emh) & exist([currdatapath,'/',filename,'.e1h'])==2)
    eval(['load ',currdatapath,'/',filename,'.e1h'])
    eval(['emh=',filename,'./40;'])
  end
  if (isempty(emv) & exist([currdatapath,'/',filename,'.e1v'])==2)
    eval(['load ',currdatapath,'/',filename,'.e1v'])
    eval(['emv=',filename,'./40;'])
  end
end

% Data loading error checks
  
if isempty(foo)
  disp('error in loading: .s file not found')
  ERROR=1;
end
if isempty(emt)
  disp('error in loading: .tim file not found')
  ERROR=1;
end      
if isempty(emh)
  disp('error in loading: .e1h file not found')
  ERROR=1;
end      
if isempty(emv)
  disp('error in loading: .e1v file not found')
  ERROR=1;
end      
		
if ERROR==0
		
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

	%%%%% Get Paradigm Number and FP Coordinates %%%%%
	for dpi=1:length(datapath)
	  fid=fopen([char(datapath(dpi)),'/',filename,'.s']);
	  if (fid > 1)
	    break;
	  end
	end
	if fid==-1
		disp('error in loading:  cannot open .s file')
	else
		index=0;fpx=[];fpy=[];
		while index==0
			line=fgets(fid);
			if line==-1;index=1;end
			if length(line)>8
				index=strcmp(line(1:8),'paradigm');
				if index==1
					paradigm=str2num(line(10:length(line)));
				end
			end
		end
		if isempty(paradigm)
			disp('error in loading:  no paradigm number found')
		else
		
			%%% Kluge for GregH's data
			if paradigm==502
				raster=raster(:,3:size(raster,2));
			end
		
			if exist(['fpstrings',num2str(paradigm)])~=1
				disp('error in loading:  default fp strings not found for this paradigm')
			else
				eval(['fpstrings=fpstrings',num2str(paradigm),';'])
				fpxstring=char(fpstrings(1));
				fpystring=char(fpstrings(2));
				index=0;
				while index==0
					line=fgets(fid);
					if line==-1;index=1;end
					if length(line)>length(fpxstring)
						index=strcmp(line(1:length(fpxstring)),fpxstring);
						if index==1
							fpx=str2num(line(length(fpxstring)+1:length(line)))/10;
							if paradigm==363;fpx=fpx*(-1);end		% kluge
							if paradigm==364;fpx=fpx*(-1);end
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
				
				if isempty(fpx) | isempty(fpy)
					disp('error in loading:  fpx or fpy never found')
				else
					
					%%%%% Normalizing Data %%%%%
					emh=emh-fpx;
					emv=emv-fpy;

				end
			end
		end
		fclose(fid);clear fid
		clear index fpx fpy fpxstring fpystring line fpstrings
	end
end	

%%%%%%%%%%%%%%%%%%%%%%%%%%%							
%%%%% Loading Indices %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ERROR==0

if ~isempty(paradigm)
	if exist([indpath,'/indices_',num2str(paradigm),'.m'])==2
		eval(['cd ',indpath])
		eval(['indices_',num2str(paradigm),';'])
		eval(['cd ',fxnpath])

		s=who('SIN_*');
		indexcols=zeros(size(eventnames,1),1);
		for x=1:size(s,1)
			eval(['sinval=',char(s(x)),';'])
			if ~isempty(sinval)
				for y=2:size(eventnames,1)
					if strcmp(char(s(x)),char(eventnames(y)))
						indexcols(y)=sinval;
					end
				end
			end
		end
		
	else
		disp('error in loading:  indices not available')
	end		
else
	disp('error in loading:  indices not available')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Calculating Data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('   calculating data')
	
pemt=emt;
pemh=emh;
pemv=emv;
msac=[];
pmsac=[];
	
trunflg=indexcols(truncatecol);
if ~isempty(algsetname)
	msacflg=1;
else
	msacflg=0;
end
	
if trunflg

	if isempty(datindex)
		disp('error in calculating:  cannot truncate because .index file empty')
	elseif size(datindex,1)~=size(pemt,1)
		disp('error in calculating:  number of trials in .index and .tim not equal')
	else
		
		%%%%% Truncating Data %%%%%
		validvalues=datindex(:,trunflg);
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
	
end
	
if regflg

	%%%%% Regularizing Data %%%%%
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
	
if msacflg

	%%%%% Getting Microsaccades %%%%%
	algparams=[];
	if ~isempty(algsetname)
		if exist([algsetpath,'.mat'])~=2
			disp('error in calculating:  algset file not found')
		else
			eval(['load ',algsetpath])
			if exist([algsetname])~=1
				disp('error in calculating:  algset variable not found')
			else
				eval(['algparams=',algsetname,';'])
			end
		end
	end
	padding=nan;
	msac=msacalg(pemt,pemh,pemv,algparams,padding);
	pmsac=msac;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Evaluating Special Functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for x=1:size(fxncalls,1)
	eval([fxncalls{x,1},'=',fxncalls{x,2},';'])
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Setting Up Outputs %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for x=1:size(outputlist,1)
	eval(['varout',num2str(x),'=',outputlist{x,1},';'])
end

	
eval(['cd ',currpath])