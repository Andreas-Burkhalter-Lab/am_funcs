function [result]=bintowave(input_file)
%takes a .bin file and makes a .wav file
fid=fopen(input_file);
if fid==-1
  error('cant open file');
end
if ~strcmp(input_file(end-3:end),'.bin')
  error('input must be a .bin file')
fid=fopen(input_file);
else
 header = ReadAIHeaderWU1(fid);
fseek(fid,header.headersize,-1);
  

y = fread(fid,inf,'uint16');
%put on a -1 1 scale with [0 65535]
%first put on a 0 2 scale
yy=(y*2/65535)-1;
fclose(fid);
wavwrite(yy,250000,16,input_file(1:end-4));
result=1;
end
%%test it
 %sngparms.plot = false;sngparms.nfreq = 256; sngparms.freqrange = [25000 110000];sngparms.threshold = 4000;
 %sound2sng('m45.wav',sngparms,'m45from.sng');
 %spsngplot('m45from.sng')
 
 