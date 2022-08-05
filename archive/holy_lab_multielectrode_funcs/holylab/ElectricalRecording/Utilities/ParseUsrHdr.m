function [txtout,vlabel] = ParseUsrHdr(usrhdr)
% ParseUsrHdr: interpret the text placed in the AI file header
% [txtout,vlabel] = ParseUsrHdr(usrhdr)
% txtout is a cell array, where txtout{i} is one return-delimited
% line of the user header
% vlabel is a cell array from 1:12, where each line contains the text
% describing the corresponding valve contents
%        (it is assumed that the valve number appears at the beginning
%         of the relevant line)
s = char(usrhdr);
% Figure out if it's a mac or a unix file
ret = sprintf('\n');
if isempty(findstr(ret,s))
  ret = sprintf('\015');  % Try a mac \r
  if isempty(findstr(ret,s))
      ret = '\n'; % Literal \n
  end
end
%[txtout{1},r] = strtok(s,ret);
%while (~isempty(r))
%        s = r;
%        [txtout{end+1},r] = strtok(s,ret);
%        if (isempty(txtout{end}))
%                txtout(end) = [];
%        end
%end
% Break into lines
retindx = [strfind(s,ret),length(s)+1];
txtout = {};
for i = 1:length(retindx)-1
    txtout{i} = s(retindx(i)+length(ret):retindx(i+1)-1);
end
%for i = 1:length(txtout)
%        disp(txtout{i});
%end
if (nargout > 1)
        vlabel = cell(1,16);
        for i = 1:length(txtout)
                vnum = sscanf(txtout{i},'%d');
                if (isnumeric(vnum) & length(vnum) == 1 & vnum > 0)
                        eindx = findstr('=',txtout{i});
                        if (length(eindx) == 1)
                                eindx = eindx+1;
                                while isspace(txtout{i}(eindx))
                                        eindx = eindx+1;
                                end
                                vlabel{vnum} = txtout{i}(eindx:end);
                        else
                                %vlabel{vnum} = txtout{i};
                        end
                end
        end
end
