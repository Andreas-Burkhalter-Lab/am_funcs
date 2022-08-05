function usave_out = fill_bad_frames(smm_in, usave_in)
% fill_bad_frames will attempt to identify and interpolate bad frames from a .imagine/.cam field
% Syntax: usave_out = fill_bad_frames(smm_in, usave_in)
%
% Options:
%  lat_interp(default [3]): limits the lateral search for bad frame interpolation.
%
% Imaging data set should be manually truncated if there are bad tail
% acquisition stacks in order to avoid crashing this utility and
% over-interpolating data.
%
% Note that method.fix_range and method.fix_index can be supplied as special cases for
% interpolation of data.  An example is when a known manual translation
% occurs during a live recording, a smear may need to be replaced with
% some adjacent stacks.  

% See also multigrid_registration_stepper

% Copyright 2010 Julian P. Meeks (Tim Holy Laboratory)

size_smm = double(smm_in.size);
smm_sz = smm_in.size;
h = smm_in.header;
padmat = ones(1, size_smm(2), size_smm(3), 'single');
if isfield(usave_in.method,'bad_stacks')
    usave_in.bad_stacks = usave_in.method.bad_stacks;
end
if isfield(usave_in.method,'bad_frames')
    usave_in.bad_frames = usave_in.method.bad_frames;
end
if isfield(usave_in.method, 'lat_interp')
  lat_interp = usave_in.method.latinterp;
else
  lat_interp = default(3);
end
% if isfield(usave_in.method, 'is_good')
%   usave_in.is_good = usave_in.method.is_good;
% end
% if isfield(usave_in.method, 'fix_range')
%   usave_in.fix_range = usave_in.method.fix_range;
% end
% if isfield(usave_in.method, 'fix_index')
%   usave_in.fix_index = usave_in.method.fix_index;
% end
usave_in = default(usave_in,'bad_stacks',[]);
usave_in = default(usave_in,'bad_frames',[]);
usave_in = default(usave_in,'bad_stack_file',[]);
usave_in.method = default(usave_in.method,'is_good', []);
usave_in.method = default(usave_in.method,'fix_range', []);
usave_in.method = default(usave_in.method,'fix_index', []);
range_stacks = [];
range_index = [];
usave_in.method = default(usave_in.method, 'thresh',(smm_sz(1)+smm_sz(2))/2);
if logical(size(usave_in.method.fix_range,2))
  for i = 1:size(usave_in.method.fix_range,2)
    c = size(range_index,2);
    range_index((c+1):(c+size(usave_in.method.fix_range{i},2))) = i;
    range_stacks = [range_stacks usave_in.method.fix_range{i}];
  end
  usave_in.bad_stacks = union(usave_in.bad_stacks, range_stacks);
end
    
while any([isempty(usave_in.bad_stacks) isempty(usave_in.bad_frames) isempty(usave_in.bad_stack_file)])
    if isempty(usave_in.bad_stacks)
        fprintf('''fill_bad_frames'' was invoked without supplying any bad_stacks or frames, will try to find them for you...\n');
        badstack = []; %badstack is used only in this first condition of the while statement
        for i = 1:smm_sz(4)
            thisstack = single(smm_in(:,:,:,i));
            fprintf('%d...',i);
            badframe = NaN;
            findzeros = find(nansum(nansum(thisstack,2),1)==0);
            if ~isempty(findzeros)
                if isempty(badstack)
                    badstack = i;
                else
                    badstack = [badstack i];
                end
            end
                    
%             for j = 1:smm_sz(3)
%                 thisframe = double(smm_in(:,:,j,i));
%                 if nansum(nansum(thisframe,2),1)==0
%                     if isempty(badstack)
%                         badstack(1) = i;
%                     else
%                         badstack = [badstack i];
%                     end
%                 end
%             end
        end
        usave_in.bad_stacks = union(unique(badstack),range_stacks);
        usave_in.bad_stacks = setdiff(usave_in.bad_stacks, usave_in.method.is_good);
    elseif isempty(usave_in.bad_frames)
        thresh = usave_in.method.thresh;
        bframe{1} = [];
        badstacks = usave_in.bad_stacks;
        n_badstacks = size(badstacks,2);
        for i = 1:n_badstacks
            badframe = NaN;
            left = 1;
            right = 1;
            % Need to identify adjacent good stacks.
            if ismember((badstacks(i)-1),badstacks)
              c = 1;
              while c<lat_interp && ismember((badstacks(i)-c),badstacks)
                if ~ismember((badstacks(i)-c-1),badstacks)
                  left = c+1;
                  c = NaN;
                else
                  c = c+1;
                end
              end
            end
            if ismember((badstacks(i)+1),badstacks)
              c = 1;
              while c<lat_interp && ismember((badstacks(i)+c),badstacks)
                if ~ismember((badstacks(i)+c+1),badstacks)
                  right = c+1;
                  c = NaN;
                else
                  c = c+1;
                end
              end
            end
            if logical(intersect(badstacks(i),range_stacks))
              badframe = [1:smm_sz(3)];            
            else
              for j = 1:smm_sz(3)
                thisframe = double(smm_in(:,:,j,badstacks(i)));
                if nansum(nansum(thisframe,2),1)==0
                    badframe(end+1) = j;
                else  
                  
                  thismin = median(min(thisframe,[],1));
                  thismean = nanmean(double(smm_in(:,:,j,[badstacks(i)-left badstacks(i)+right])),3);  % recently changed, may need to adjust as was along 4th dimension
                  thisstd = nanstd(double(smm_in(:,:,j,[badstacks(i)-left badstacks(i)+right])),0,4);
                  theseoff = (thisframe.^2)>(thismean+thismin).^2;
                  thesetoo = ((thismean-thismin).^2)>(thisframe.^2);
                  theseoff = theseoff + thesetoo;
                  if sum(sum(theseoff,2),1) > thresh
                      badframe(end+1) = j;
                  end
                end
              end
            end            
            badframe(isnan(badframe)) = [];
            if ~isempty(badframe)
                if isempty(bframe{1})
                    bframe{1} = badframe;
                else
                    bframe{end+1} = badframe;
                end
            else
                usave_in.bad_stacks(i) = NaN;
            end
        end
        usave_in.bad_stacks(isnan(usave_in.bad_stacks))=[];
        usave_in.bad_frames = bframe;
    elseif isempty(usave_in.bad_stack_file)
        badstacks = usave_in.bad_stacks;
        n_badstacks = size(badstacks,2);
        temp = strfind(usave_in.infile, '.imagine');
        if ~isempty(usave_in.basefilename)
          basefilename = usave_in.basefilename;
        elseif ~isempty(temp)
            basefilename = usave_in.infile(1:temp-1);
        elseif ~isempty(strfind(usave_in.infile, '.cam'))
            temp = strfind(usave_in.infile, '.cam');
            if ~isempty(temp)
                basefilename = usave_in.infile(1:temp-1);
            end
        else
            basefilename = usave_in.infile;
        end
        usave_in.bad_stack_file = [basefilename '_bad_stack.cam'];
        [fid,msg] = fopen(usave_in.bad_stack_file, 'w');
        if (fid < 0)
            fprintf('%s\n',msg);
            error(['Can''t open file ' savedir savefile ' for writing. Do you have permission?']);
        end
        for i = 1:n_badstacks
            thisstack = double(smm_in(:,:,:,badstacks(i)));
            badframe = usave_in.bad_frames{i};
            %  We want to check adjacent frames and assign first good
            %  frames to be used.
            left = 1;
            right = 1;
            if ismember((badstacks(i)-1),badstacks)
              c = 1;
              while c<lat_interp && ismember((badstacks(i)-c),badstacks)
                if ~ismember((badstacks(i)-c-1),badstacks)
                  left = c+1;
                  c = NaN;
                else
                  c = c+1;
                end
              end
            end
            if ismember((badstacks(i)+1),badstacks)
              c = 1;
              while c<lat_interp && ismember((badstacks(i)+c),badstacks)
                if ~ismember((badstacks(i)+c+1),badstacks)
                  right = c+1;
                  c = NaN;
                else
                  c = c+1;
                end
              end
            end
            if logical(intersect(badstacks(i),range_stacks))              
              thisstack(:,:,:) = nanmean(double(smm_in(:,:,:,[usave_in.method.fix_index{range_index(range_stacks==badstacks(i))}(1) ...
                usave_in.method.fix_index{range_index(range_stacks==badstacks(i))}(end)])),4);
            else for j = 1:size(badframe,2)
                thisstack(:,:,badframe(j)) = nanmean(double(smm_in(:,:,badframe(j),[badstacks(i)-left badstacks(i)+right])),4);
              end
            end
            if strcmp(h.camera, 'DV8285_BV'), 
              thisstack = cat(1,thisstack,padmat);
              count = fwrite(fid, uint16(thisstack), 'uint16');              
            else
              count = fwrite(fid, uint16(thisstack), 'uint16');
            end
            if (count < prod(size_smm(1:end-1)))
                error('Wrote fewer pixels than expected, quitting');
            end
        end
    end
end
usave_out = usave_in;
fclose(fid);
% TODO: Display a popup demanding the user hit "Okay" to set up the
% .imagine file that corresponds to the # of bad_stacks (or just do it with
% a script to replace the text of a new imagine file automatically??
end