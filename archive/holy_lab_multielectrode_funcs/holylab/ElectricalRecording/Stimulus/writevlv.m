function writevlv(fileout,stimfiles,stim)
% writevlv: write the time course of the stimulus in the format of CalcStim
% writevlv(outfilename,files,stim)
% where
%        outfilename is the base name for output
%        files is a cell array of filenames
%        stim is a cell array of 2-by-n matrices, containing the
%                valve number and scan number of each transition
% 

% Used to write both .sat (= Stimulus Artifact Time) and .vlv files,
% but now there are no stimulus artifacts so the .sat file is not written
[pathstr,basename,ext] = fileparts(fileout);
if(isempty(pathstr))
  pathstr=pwd; % @jason: in case user gives a name like 'stim.vlv' w/o path
end

if strcmp(ext,'.vlv')
  fileout = [pathstr filesep basename];
end
[fid,message] = fopen([fileout,'.vlv'],'wt');
if (fid < 0)
        error(message);
end
for i = 1:length(stimfiles)
        fprintf(fid,'%s {',stimfiles{i});
        fprintf(fid,'  %d %d',round([stim{i}(1,:);stim{i}(2,:)])); % @jason: "%d of double" outputs in double format sometimes
        fprintf(fid,'}\n');
end
fclose(fid);
return
[fid,message] = fopen([fileout,'.sat'],'wt');
if (fid < 0)
        error(message);
end
for i = 1:length(stimfiles)
        fprintf(fid,'%s {',stimfiles{i});
        t = stim{i}(2,:);
        killi = find(t(1:end-1)+600 > t(2:end));
        t(killi) = [];
        fprintf(fid,'  %d %d',[t(1:end-1)+600;t(2:end)]);
        fprintf(fid,'}\n');
end
fclose(fid);
