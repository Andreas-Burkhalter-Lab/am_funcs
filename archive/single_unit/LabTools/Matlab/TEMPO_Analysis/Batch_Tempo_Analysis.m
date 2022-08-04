function Batch_Tempo_Analysis

TEMPO_Defs;
ProtocolDefs;
Path_Defs;

batchfiledir = 'Z:\Data\Tempo\Batch Files\'
filename = input('Enter Batch File Name: ', 's');
close_flag = input('Do you wish to close all windows (0=N, 1=Y)? ');
print_flag = input('Do you wish to print figures (0=N, 1=Y)? ');

status = ['Loading batch file: ' batchfiledir filename];
disp(status);

filename = [batchfiledir filename];
fid = fopen(filename);

line = fgetl(fid);
while (line ~= -1)
   
   if (line(1) ~= '%')
	   spaces = isspace(line);
		space_index = find(spaces);

		%get path / file
		PATH = line(1:space_index(1) - 1);
		FILE = line(space_index(1) + 1:space_index(2) - 1)

		%get analysis type
		%special commands for analysis read
		i = 2;
		TempAnalysis = '';
		while(line(space_index(i)-1) ~= '''')
		      i = i+1;
		end
   
		TempAnalysis = strcat(TempAnalysis, line(space_index(2) + 2:space_index(i)-2));
		Analysis = {};
		Analysis{1} = TempAnalysis;

		%get beginning trial
		BegTrial = line(space_index(i) + 1:space_index(i+1) - 1);
		BegTrial = str2num(BegTrial);

		%get ending trial
		EndTrial = line(space_index(i+1) + 1:space_index(i+2) - 1);
		EndTrial = str2num(EndTrial);

		%get startcode
		StartCode = line(space_index(i+2) + 1:space_index(i+3) - 1);
		StartCode = str2num(StartCode);

		%get startoffset
		StartOffset = line(space_index(i+3) + 1:space_index(i+4) - 1);
		StartOffset = str2num(StartOffset);

		%get stopcode
		StopCode = line(space_index(i+4) + 1:space_index(i+5) - 1);
		StopCode = str2num(StopCode);

		%get stopoffset
		StopOffset = line(space_index(i+5) + 1:space_index(i+6) - 1);
		StopOffset = str2num(StopOffset);

		%get select_data - good trials or bad trials
		good_flag = line(space_index(i+6) + 1:space_index(i+7) - 1);
		good_flag = str2num(good_flag);

		%get spikechannel
		SpikeChan = line(space_index(i+7) + 1:length(line));
		SpikeChan = str2num(SpikeChan);

		%load data file
		[return_val, g_data, b_data] = LoadTEMPOData(PATH,FILE);
		if (good_flag)
			select_data = g_data;
		else
			select_data = b_data;
	   end
   	if BegTrial == -1
      	BegTrial = 1;		
	   end
   
   	if EndTrial == -1
      	EndTrial = size(g_data.event_data,3);		
	   end
   
		%get protocol type
		Protocol = g_data.one_time_params(PROTOCOL);		%get the protocol string descriptio

		%analyze data
		batch_flag = 1;
	   Analysis_Switchyard(select_data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);
   
	   if print_flag
   	   print;
	   end
   
   	if close_flag
	      close all;
      end
   end % if (line(1) ~=...
  	line = fgetl(fid);
   
end %while...

fclose(fid);