function SaveLB(filein,fileout)
% SaveLB: Write behaviors logged by LogBehavior to a file
% SaveLB(filein,fileout), where filein is the name of the .bin file
% and fileout is the name of the text file for output.

v = LogBehavior(filein);
% See if this file exists already
[fido,message] = fopen(fileout,'rt');
newfile = 0;
if (fido < 0)
        newfile = 1;
else
        fclose(fido);
end
[fidw,message] = fopen(fileout,'at');
if (fidw < 0)
        error(message);
end
if (newfile)        % If new file, prepend names of behaviors
        behavs = LogBehavior;
        fprintf(fidw,'%d\n',length(behavs));
        for i = 1:length(behavs)
                fprintf(fidw,'%s\n',behavs{i});
        end
end
fprintf(fidw,'%s\n',filein);
for i = 1:length(v)
        fprintf(fidw,'%d\n',size(v{i},2));
        fprintf(fidw,'%3.2f  %d\n',v{i});
end
fclose(fidw);
