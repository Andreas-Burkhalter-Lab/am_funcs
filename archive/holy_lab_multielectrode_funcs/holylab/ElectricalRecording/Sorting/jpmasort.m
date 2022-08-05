function asort(outname,innamebase,op)

% Sets up and runs the autosort so that all you have left
% to do is run "cass" and select whatever filename you gave it
% for output (eg file_out), then run "cass_sort_apply('file_out')"
%
% if only one input is given, all .ssnp files in the directory are 
% used, and the one input is used as the outname; if "clean" is 
% given as the inname, the cleaned versions of all .ssnp files 
% in the directory are used
%
% if a cell array of file names are given instead of "innamebase" 
% it will just use those names

if nargin < 2
  infiles = dirbyname('*.ssnp');
elseif strcmp(innamebase,'clean')
  infiles = dirbyname('*clean.ssnp');
elseif strcmp(innamebase,'cleanECG')
  infiles = dirbyname('*clean.ECG_ssnp');
elseif strcmp(innamebase,'ECG')
    infiles = dirbyname('*.ECG_ssnp');
elseif iscell(innamebase)
  infiles = innamebase;
else
  infiles = dirbyname([innamebase '*.ssnp']);
end

op.n_replicates = 1;
%op.Rfactor = [1];
%op.Rfactor = [.8 .95 1.15 1.5 2 3];
op = default(op,'max_snip_memsize',5e7);
op = default(op,'max_snips_to_cluster',Inf);
op = default(op,'t2Vfunc',@normalized_t2V);
op = default(op,'cluster_func',@msams);
op = default(op,'projection_func',@pca_sort_fixed);
op = default(op,'interactive_kmeans',1);
% If they were set the opt passed to asort, 't2V_tmult' and 'n_landmarks'
% will also be passed along
op = default(op,'n_landmarks',20000);
op = default(op,'factor',3);     % for msams_step1
op = default(op,'min_to_check',3);% "
% to increase # of principal components (JPM 2007_02_28)
op = default(op,'pca_n_components', 3);

sh = snipfile2sortheader(infiles);
% Please run with the following two enabled for a while; it will help
% shake out rare bugs
dbstop if error
dbstop if warning 'MATLAB:divideByZero'
% Now do the clustering
autosort(sh,outname,op)



