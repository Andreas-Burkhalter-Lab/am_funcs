function result=sortfilebytime(files)
% sort files by modification time (from the earliest to the latest)
% result=SORTFILEBYTIME(files)
% pre: 
%    files: a cell array of filenames
% post:
%    result: sorted file names
   
   % get the date for each file
   for idx=1:length(files)
      dd(idx)=dir(files{idx});
      tt(idx)=datenum(dd(idx).date);
   end
   
   [foo, indices]=sort(tt);
   result=cell(1, length(files));
   [result{:}]=deal(files{indices}); % see help on LISTS for meaning of [aCellArray{:}]
   % result{:} % print out the result