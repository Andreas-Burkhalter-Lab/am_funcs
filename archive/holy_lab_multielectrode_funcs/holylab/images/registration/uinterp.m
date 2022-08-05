function usave_out = uinterp(smm_in, usave_out, options)
% uinterp interpolates "u" values from image registration

% Copyright 2010 Julian P. Meeks (Tim Holy Laboratory)

%% set options
options = default(options,'startstack',1);
options = default(options,'mode','linear');
usave_out.method = default(usave_out.method,'save_new_u',true);
options.save_new_u = usave_out.method.save_new_u;
usave_out.method = default(usave_out.method,'output','none');

%% identify limits

startstack = options.startstack;
laststack = smm_in.size; laststack = laststack(4);
u_dims = size(usave_out.u_matrix{1}); u_dims = u_dims(1:3);
[sortstacknum sortstack_idx] = sort(usave_out.stacknum);

writecam = false;
if ~isempty(strmatch('none',usave_out.method.output))
    writecam = false;
elseif ~ischar(usave_out.method.output)
    error('''method.output'' must contain a string variable.\n');
else 
    temp = findstr('.imagine',usave_out.infile);
    if ~isempty(temp)
        base = usave_out.infile(1:temp(end)-1);
        outfile = [base '_' usave_out.method.output '.cam'];
    end
    temp = findstr('.cam', usave_out.infile);
    if ~isempty(temp)
        base = usave_out.infile(1:temp(end)-1);
        outfile = [base '_' usave_out.method.output '.cam'];
    end
    if ~exist('base', 'var')
        base = usave_out.infile;
        outfile = [base '_' usave_out.method.output '.cam'];
    end
    writecam = true;
end


% gather stack info
stack_size = smm_in.size;
header_smm = smm_in.header;
stimuli = header_smm.stim_lookup;
Number_of_stacks = stack_size(end);
stimuli = stimuli(1:Number_of_stacks);


% open the new .cam file
if writecam == true
    fid = fopen(outfile,'w');
    fprintf('Writing stacks:\n')
end


%% apply interpolation
% check to see if we're going to save the uvalues (only okay at coarser grids)
if options.save_new_u
    utempmatrix = cell([1 laststack+1-startstack]);
end
count = 1;
if ~isempty(strmatch('linear', options.mode, 'exact'))
    thisstart = startstack;
    for i = 1:size(usave_out.stacknum,2)+1
        if i == 1
            interpn = sortstacknum(i)-thisstart;
            divn = sortstacknum(i+1)-sortstacknum(i);
            gradient = (usave_out.u_matrix{sortstack_idx(i)}-usave_out.u_matrix{sortstack_idx(i+1)})/divn;
            for j = 1:interpn
               utempmatrix{count} = usave_out.u_matrix{sortstack_idx(i)}-gradient*(interpn+1-j);
               if writecam
                   warpnwrite(fid, smm_in, utempmatrix{count}, sortstacknum(i)-(interpn+1-j), usave_out)
               end
               if options.save_new_u
                   count = count+1;
               end
            end
        elseif i > size(usave_out.stacknum,2)
            interpn = laststack - sortstacknum(i-1);
            divn = sortstacknum(i-1)-sortstacknum(i-2);
            gradient = (usave_out.u_matrix{sortstack_idx(i-1)}-usave_out.u_matrix{sortstack_idx(i-2)})/divn;
            for j = 1:interpn+1
               utempmatrix{count} = usave_out.u_matrix{sortstack_idx(i-1)}+gradient*(j-1);
               if writecam
                   warpnwrite(fid, smm_in, utempmatrix{count}, sortstacknum(i-1)+j-1, usave_out)
               end
               if options.save_new_u
                   count = count+1;
               end
            end
        else
            interpn = sortstacknum(i)-sortstacknum(i-1);
            if isempty(usave_out.u_matrix{sortstack_idx(i)})
                temp = zeros(size(usave_out.u_matrix{i-1}));
                gradient = (temp-usave_out.u_matrix{sortstack_idx(i-1)})/interpn;
            elseif isempty(usave_out.u_matrix{sortstack_idx(i-1)})
                temp = zeros(size(usave_out.u_matrix{i}));
                gradient = (usave_out.u_matrix{sortstack_idx(i)}-temp)/interpn;
            else
                gradient = (usave_out.u_matrix{sortstack_idx(i)}-usave_out.u_matrix{sortstack_idx(i-1)})/interpn;
            end
            for j = 1:interpn
               if isempty(usave_out.u_matrix{sortstack_idx(i-1)})
                   utempmatrix{count} = temp+gradient*(j-1);
               else
                   utempmatrix{count} = usave_out.u_matrix{sortstack_idx(i-1)}+gradient*(j-1);
               end
               if writecam
                   warpnwrite(fid, smm_in, utempmatrix{count}, sortstacknum(i-1)+j-1, usave_out)
               end
               if options.save_new_u
                   count = count+1;
               end
            end
        end
    end
elseif ~isempty(strmatch('spline', options.func, 'exact'))
% NOT YET IMPLEMENTED
    error('not yet implemented.\n');
else
    error('not yet implemented.\n');
end

%% set output
if isfield(usave_out,'u_err')
    usave_out.u_err = [];
end
if isfield(usave_out, 'orig_err');
    usave_out.orig_err = [];
end
if isfield(usave_out,'n_cycles')
    usave_out.n_cycles = [];
end
if options.save_new_u
    usave_out.u_matrix = utempmatrix;
else
    usave_out.u_matrix = [];
    usave_out.u_matrix = 'linear_uinterp_not_saved';
end
clear utempmatrix;
usave_out.stacknum = 1:size(usave_out.u_matrix,2);
usave_out.method.mode = options.mode;
usave_out.method.save_new_u = options.save_new_u;
%% ------------------------
end