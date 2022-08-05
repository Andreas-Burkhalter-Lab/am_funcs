 
 fid = fopen('/home/eric/t.bin', 'r');
 [header] = ReadAIHeader(fid);
 fclose(fid);
%  
%  fid = fopen('/home/eric/d.bin', 'r');
%  header = ReadAIHeader(fid);
%  factor(header.nscans/2*header.nfreq)

%  fid = fopen('/home/eric/m1sF.det', 'r');
%  ReadSensor(fid);
%  fclose(fid);

%  fid = fopen('/usr/lab/matlabfunc/Whistles/E8_5.sng', 'r');
% % timeRange = [0 7];
%  header = ReadSonogram(fid);
 
