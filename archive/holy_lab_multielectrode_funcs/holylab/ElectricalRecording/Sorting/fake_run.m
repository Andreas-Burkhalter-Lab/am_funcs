op = struct('n_replicates',2,'Rfactor',[1 1.5 2 3 4 5 6 8 10]);
fls = dirbyname('fake*.ssnp');
sh = snipfile2sortheader(fls);
% Please run with the following two enabled for a while; it will help
% shake out rare bugs
dbstop if error
dbstop if warning 'MATLAB:divideByZero'
% Now do the clustering
autosort(sh,'sort1',op)

% When it's finished: type 'cass' and choose the sort1 directory

