function [stimname,stimconc] = parsestimnames(valvelabels)
% PARSESTIMNAMES: return identity & concentration from valvelabel
% Syntax:
%  [stimname,stimconc] = parsestimnames(valvelabels)
% The valvelabel strings are parsed for '/' as a separator between the name
% and concentration.
% For example, if
%   valvelabels = {'F/100','F/1000'}
% then a call to PARSESTIMNAMES yields
%   stimname = {'F','F'}
%   stimconc = [0.01 0.001]

for i = 1:length(valvelabels)
    % First get rid of spurious #=KCl cases
    vl = valvelabels{i};
    eqindx = findstr(vl,'=');
    if ~isempty(eqindx)
        tok = [num2str(sscanf(vl,'%d')) '='];
        if strncmp(tok,vl,length(tok))
            vl = vl(length(tok)+1:end);
        end
    end
    [stimname{i},rest] = strtok(vl,'/');
    tmp = str2num(rest(2:end));
    if ~strcmp([stimname{i} '/' num2str(tmp)],vl)
        error('Can''t parse all names properly, some don''t conform to pattern');
    end
    stimconc(i) = 1/tmp;
end
    