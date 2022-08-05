
% this script registers all the stacks

clear all

profile on
datadir = '/sataext3/dturaga/2006-02-02';
savedir = '/sataext3/dturaga/registered';
datafile = '2006-02-02-2';
savefile = [datafile '_reg.cam'];
savefile_g = [datafile '_g.cam'];
base_stack_index = 1;

% Set the decimation schedule for the images that will be used to estimate
% the warp
base_decimation_schedule = [2 2 1; 2 2 1];

% Create the decimation schedule needed to get the warp grid size from the
% (reduced) images.
options = struct;
options.g_decimate = [3 3 1; 2 2 2; 2 2 2; 2 2 2];
options.step2_thresh = 1e-6;  % convergence criterion
  
cd /sataext3/dturaga/2006-02-02
smm = stackmm(datafile);
stack_size = smm.size;
Number_of_stacks = stack_size(end);

if exist([savedir filesep savefile],'file')
    response = input('File exists. Replace or append? (r/a) [a]:', 's');
    if strcmp(lower(response),'r')
        delete([savedir filesep savefile]);
    end
end

if exist([savedir filesep savefile_g],'file')
    response = input('File exists with g values. Replace or append? (r/a) [a]:', 's');
    if strcmp(lower(response),'r')
        delete([savedir filesep savefile_g]);
    end
end

[fid,msg] = fopen([savedir filesep savefile],'a');
if (fid < 0)
    fprintf('%s\n',msg);
    error(['Can''t open file ' savedir filesep savefile ' for writing. Do you have permission?']);
end

[fid_g,msg] = fopen([savedir filesep savefile_g],'a');
if (fid < 0)
    fprintf('%s\n',msg);
    error(['Can''t open file ' savedir filesep savefile_g ' for writing. Do you have permission?']);
end

stack_base = smm(:,:,:,base_stack_index);
stack_base = single(stack_base);
for decIndex = 1:size(base_decimation_schedule,1)
    stack_base = imreduce(stack_base, base_decimation_schedule(decIndex,:));
end
psi1 = sqrt(stack_base); clear stack_base
options.output_size = size(psi1);  % size to expand g to

n_decimation_steps = size(options.g_decimate,1);
% Now actually run this decimation sequence to see how big g needs to be
psi1tmp = psi1;
for decIndex = 1:n_decimation_steps
  psi1tmp = imreduce(psi1tmp,options.g_decimate(decIndex,:));
end
g_sz = size(psi1tmp);
g = register_g0(g_sz)   % Let the user see the output

fprintf('Progress (out of %d stacks):\n',Number_of_stacks);
mu = 0.1;
err = zeros(1,Number_of_stacks);

for i = base_stack_index+1:Number_of_stacks
  fprintf('..%d..',i);
    
    stack_temp = smm(:,:,:,i);
    stack_temp = single(stack_temp);
    for decIndex = 1:size(base_decimation_schedule,1)
        stack_temp = imreduce(stack_temp, base_decimation_schedule(decIndex,:));
    end
    psi2 = sqrt(stack_temp); clear stack_temp
    mu = max(mu,0.01);
    [g,psig,mu,err(i)] = register_nonrigid(psi1,psi2,g,0,mu,options);

% recreating the stacks because the original ones were image reduced
    clear psig
    stack_temp = smm(:,:,:, i);
    stack_temp = single(stack_temp);
    
    g_hires = register_expandg(g,size(stack_temp),'iminterp');
   
    img = register_warp(stack_temp,g_hires);

    size_smm = double(smm.size);
    padmat = ones(1, size_smm(2), size_smm(3), 'single');
    img = cat(1, img, padmat);
    
    % Save the result
    count = fwrite(fid,uint16(img),'uint16');
    if (count < prod(stack_size(1:end-1)))
        error('Wrote fewer pixels than expected, quitting');
    end
    
    g_mat = cell2mat(g);
    fwrite(fid_g, uint16(g_mat), 'uint16');
    
    clear g_hires img stack_temp
    
end



