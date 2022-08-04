function[]=compilemsac(filelist,savefilename);

% function[]=compilemsac(filelist,savefilename);
% Function (really a script) to loop through data getting msac
%	stats for population curves.

savedir='E:/data/matlab';

doneflag=0;

saveevents=[{'SIN_FPON'};{'SIN_FIX'};{'SIN_STON'};{'SIN_STOFF'}];
outputlist=[{'msac'};{'pemt'};{'pemh'};{'pemv'};{'datindex'};{'indexcols'};{'eventnames'}];
filenumtot=size(filelist,1);

exptflag=0;
exptstatsall=[];
trlstatsall=[];
msacstatsall=[];
intstatsall=[];

eventindex=[];

if exist([savedir,'/',savefilename,'.mat'])~=2
	eval(['save ',savedir,'/',savefilename,' exptflag exptstatsall trlstatsall msacstatsall intstatsall saveevents'])
end

while ~doneflag

	eval(['load ',savedir,'/',savefilename])
	filenum=exptflag+1;
	
	if filenum<=filenumtot

		exptflag=filenum;
		disp(['file ',num2str(filenum),' of ',num2str(filenumtot)])
		
		filename=filelist(filenum,:);
		[msac,pemt,pemh,pemv,datindex,indexcols,eventnames]=getem(filename,outputlist);
		
		if ~isempty(pemt)

			[exptstats,trlstats,msacstats,intstats]=getmsacstats(msac,pemt,pemh,pemv);
		 
		 	% Note:  the variables saved to a file are identical to the
		 	% variables that getmsacstats spits out with the file number
		 	% (row number of file name in inputted "filelist") appended
		 	% to the end of each row.  Trlstats also includes a few
		 	% columns of event times.
		 	
			if isempty(eventindex)
				for x=1:size(saveevents,1)
					eventindex=[eventindex strmatch(saveevents{x},eventnames)];
				end
			end

			eventtimes=datindex(:,indexcols(eventindex));

			exptstats=[exptstats ones(size(exptstats,1),1)*filenum];
			trlstats=[trlstats eventtimes ones(size(trlstats,1),1)*filenum];
			msacstats=[msacstats ones(size(msacstats,1),1)*filenum];
			intstats=[intstats ones(size(intstats,1),1)*filenum];
		
			if isempty(exptstats)
		
				exptstatsall=exptstats;
				trlstatsall=trlstats;
				msacstatsall=msacstats;
				intstatsall=intstats;
			
			else

				exptstatsall=[exptstatsall;exptstats];
				trlstatsall=[trlstatsall;trlstats];
				msacstatsall=[msacstatsall;msacstats];
				intstatsall=[intstatsall;intstats];
					
			end
		end

		eval(['save ',savedir,'/',savefilename,' exptflag exptstatsall trlstatsall msacstatsall intstatsall saveevents'])

	else
	
		doneflag=1;
		disp(' -- finito --')
		
	end

end

