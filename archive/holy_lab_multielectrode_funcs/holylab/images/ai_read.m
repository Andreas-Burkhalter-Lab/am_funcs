function ai_out = ai_read(filename)
% this script helps in reading .ai files
% INPUT
% filename is the name of imagine file WITHOUT .imagine extension
% OUTPUT
% ai_out is matrix, each row containing values of one analog input
% IMPORTANT, fread as 'int16'
% Diwakar Turaga, 8-21-2008
%
% see also: STACKMM

smm = stackmm(filename);
h = smm.header;
ai_num = length(h.ai_label_list);

filename_ai = [filename '.ai'];
fid = fopen(filename_ai,'r');

ai_tot = fread(fid, 'int16');
sz = length(ai_tot);

ai_out = reshape(ai_tot, [ai_num sz/ai_num]);

fclose(fid);

end





