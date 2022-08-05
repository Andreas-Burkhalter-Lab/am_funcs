function ephysout = ephysjoincells(ephysin,cellnumsin,cellnumsout)
% EPHYSJOINCELLS: join two or more cells into one
% Syntax:
%   ephysjoincells(ephysin,cellnumsin,cellnumsout)
% where
%   ephysin is an input ephys structure array;
%   cellnumsin is an input cell array, where each element of the cell array
%     specifies the list of input cell numbers combined into the
%     corresponding output cell number;
%   cellnumsout is the corresponding output cell number;
%   
% and
%   ephysout is the output ephys structure array.
%
% Example:
% Suppose a recording gives you 6 different cells (cell numbers 1-6), and
% you want to explore what happens when you join some of them:
%    ep2 = ephysjoincells(ep1,{[2 4],[3 5 6],1},[35 36 34])
% will result in the following mapping:
%    cells 2 and 4 in ep1 will be joined into a new cell, called cell 35
%      in ep2;
%    cells 3, 5, and 6 in ep2 will be joined to make cell 36 in ep2;
%    cell 1 in ep1 will be called cell 34 in ep2.
%
% See also: EPHYS.
  
% Copyright 2004 Timothy E. Holy
  
% Right now this assumes that the cell numbers are constant across all
% elements of the array, but this could clearly be changed if
% there was a need.
  ephysout = ephysin;
  % Organize information about cell numbers
  changecells = cat(2,cellnumsin{:});  % The cells whos identities will change
  changecellsindx = findainb(changecells,ephysin(1).cellnums); % Where
                                                               % these
                                                               % are in
                                                               % the list
  newcellnums = [setdiff(ephysin(1).cellnums,changecells),cellnumsout];
  newcellnumssort = sort(newcellnums);
  newcellindx = findainb(cellnumsout,newcellnumssort);  % Where to put
                                                        % the new cells
  for i = 1:length(ephysout)
    ephysout(i).cellnums = newcellnumssort;
    ct = ephysin(i).celltimes;                    % Convenient handle to data
    ephysout(i).celltimes(changecellsindx) = [];  % delete the replaced cells
    for j = 1:length(cellnumsin)
      flagempty = zeros(1,length(cellnumsin{j})); % To avoid problems
                                                  % with empty entries
      for k = 1:length(cellnumsin{j})
        flagempty(k) = isempty(ct{cellnumsin{j}(k)});
      end
      % Concatentate data from joined cells
      ephysout(i).celltimes{newcellindx(j)} = ...
          sort(cat(1,ct{cellnumsin{j}(find(flagempty == 0))}));
    end
  end
