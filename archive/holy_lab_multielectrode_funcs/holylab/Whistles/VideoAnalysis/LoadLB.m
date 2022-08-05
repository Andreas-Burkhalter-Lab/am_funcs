function [bnames,fnames,v] = LoadLB(filename)
[fid,message] = fopen(filename);
if (fid == -1)
        error(message);
end
nbehavs = fscanf(fid,'%d',1);
i = 1;
while (i <= nbehavs)
        temp = fgetl(fid);
        if (~isempty(temp))
                bnames{i} = temp;
                i = i+1;
        end
end
i = 1;
temp = fgetl(fid);
while (temp ~= -1)
        fnames{i} = temp;
        for j = 1:nbehavs
                ntrans = fscanf(fid,'%d\n',1);
                v{i,j} = fscanf(fid,'%f %d\n',[nbehavs,ntrans])';
        end
        i = i+1;
        temp = fgetl(fid);
end
fclose(fid);
