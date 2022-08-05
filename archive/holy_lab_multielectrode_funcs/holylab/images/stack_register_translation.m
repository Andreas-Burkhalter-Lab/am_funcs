function [dx_list outstack] = stack_register_translation(instack, varargin)
% STACK_RESGISTER_TRANSLATION translation-registers a stackmm object with given options
% 
% Takes input variable "instack" as a stackmm object
% Returns outputs:
%    "dx_list": a cell array of 3-vectors with same # of stacks as input
%    "outstack": registered stack with size <= size of input stack
%        (may be specified within options.mask as described below)
% VARARGIN allows a single input structure "options" with the following subfields:
%    reference_stack = integer corresponding to the stack upon which the
%      registration will be centered.  (default = 2)
%    savefile = string containing the name of a ".cam" file to save.  If
%      file exists, a prompt will ask for permission to overwrite (append
%      otherwise).  (default = no saved file)
%    mask is a structure with the following fields
%      mask.x = [xstart xend] which gives the range of accepted x-values
%      mask.y = [ystart yend] which gives the range of accepted y-values
%      mask.z = [zstart zend] which gives the range of accepted z-values
%               (frames within the stack)
%      mask.t = [tstart tend] which gives the range of accepted time values
%               (as specified by stack number, NOT minutes or seconds)
%    OPTIONS for REGISTER_TRANSLATION_MULTIRES will be passed on 
%    dx_start = same as dx_start option for register_translation_multires
%    pixel_spacing = same as pixel_spacing for register_translation_multires
%    min_pixels = same as min_pixels for register_translation_multires
%    n_vcycles = same as n_vcycles for register_translation_multires
%    display = same as display for register_translation_multires
%
% See also REGISTER_TRANSLATION_MULTIRES, STACKMM

% Copyright 2008 Julian P. Meeks (Timothy Holy Laboratory)
% Revision History:
% 2008_03_23: Wrote it (adapted from original registration_script by DT)

%% VARARGIN CHECKING:
if nargin > 1
    options = varargin{1};
    if ~isstruct(options)
        errormsg = '''options'' argument must be a structure variable.';
        error(errormsg);
    end
end

% Set defaults for stack-related options structure
ops = default(options,'reference_stack', 2);
ops = default(options,'savefile', []);
ops = default(options,'mask',[]);
% Set defaults for register_translation_multires
t_ops = default(options,'dx_start',[0 0 0]);
t_ops = default(options,'pixel_spacing',[1 1 1]);
t_ops = default(options, 'display', false);

%% Set up/check savefiles
if ~isempty(ops.savefile)
    if exist(ops.savefile, 'file')
        overwrite = questdlg([ops.savefile ' exists. Overwrite?'], 'No');
        if ~isempty(strmatch(overwrite, 'Yes'))
            delete(ops.savefile);
        elseif ~isempty(strmatch(overwrite, 'No'))
            while exist(ops.savefile, 'file')
                newfile = inputdlg('Please specify a different savefile name',...
                                        'Change file name',...
                                        1, {ops.savefile});
                ops.savefile = newfile{1};
            end
        elseif ~isempty(strmatch(overwrite, 'Cancel'))
            return;
        end
                            
    end
    fid = fopen(ops.savefile, 'a');
    
    % if ops.savefile is given as full path, strip that path
    if ops.savefile(1) == filesep && length(strfind(ops.savefile,filesep))>1
        base_path = strip_filename(ops.savefile);
    elseif ~isempty(strfind(ops.savefile,filesep))
        base_path = strip_filename(ops.savefile);
    else
        base_path = [pwd filesep];
    end
end
% ASSERT: fid must be open for writing!

% set up temporary file directory
home = pwd;
cd(base_path);
mkdir('temp');
temp_path = [base_path 'temp/'];
cd(home);

%% Set up boundaries (x, y, z, t)
stack_size = instack.size;
n_stacks = stack_size(end);
stack_header = instack.header;
stimuli = stack_header.stim_lookup;
padmat = ones(1, stack_size(2), stack_size(3), 'single');

% Check for input masks
   if ~isempty(ops.mask)&&isstruct(ops.mask)
%  X-dimension
       if isfield(ops.mask, 'x')
           if length(ops.mask.x) == 2
               if ops.mask.x(1) < 1 || ops.mask.x(2)>stack_size(1)
                   errormsg = 'options.mask.x dimension values are not valid.';
                   error(errormsg);
               end
           else
               errormsg = 'options.mask.x must contain a 2-vector.';
               error(errormsg);
           end
       else
           ops.mask.x = [1 stack_size(1)];
       end
%  Y-dimension
       if isfield(ops.mask, 'y')
           if length(ops.mask.y) == 2
               if ops.mask.y(1) < 1 || ops.mask.y(2)>stack_size(2)
                   errormsg = 'options.mask.y dimension values are not valid.';
                   error(errormsg);
               end
           else
               errormsg = 'options.mask.y must contain a 2-vector.';
               error(errormsg);
           end
       else
           ops.mask.y = [1 stack_size(2)];
       end
%  Z-dimension
       if isfield(ops.mask, 'z')
           if length(ops.mask.z) == 2
               if ops.mask.z(1) < 1 || ops.mask.z(2)>stack_size(3)
                   errormsg = 'options.mask.z dimension values are not valid.';
                   error(errormsg);
               end
           else
               errormsg = 'options.mask.z must contain a 2-vector.';
               error(errormsg);
           end
       else
           ops.mask.z = [1 stack_size(3)];
       end
% T-dimension
       if isfield(ops.mask, 't')
           if length(ops.mask.t) == 2
               if ops.mask.t(1) < 1 || ops.mask.t(2)>stack_size(4)
                   errormsg = 'options.mask.t dimension values are not valid.';
                   error(errormsg);
               end
           else
               errormsg = 'options.mask.t must contain a 2-vector.';
               error(errormsg);
           end
       else
           ops.mask.t = [1 stack_size(4)];
       end
   else
       ops.mask.x = [1 stack_size(1)];
       ops.mask.y = [1 stack_size(2)];
       ops.mask.z = [1 stack_size(3)];
       ops.mask.t = [1 stack_size(4)];
   end
% ASSERT: if any mask options supplied, they are maintained.  If any are
% supplied, the remaining options are given defaults if not specified
% explicitly in the options.mask structure
   
% Check reference stack (and that it isn't in-line to be masked out!
  if ~isempty(ops.reference_stack)
      if ops.reference_stack < 2 || ops.reference_stack > stack_size(4)
          errormsg = ['options.reference stack is out-of bounds, choose a value 2-' num2str(stack_size(4))];
          error(errormsg);
      elseif ops.reference_stack < ops.mask.t(1) || ops.reference_stack > ops.mask.t(2)
          errormsg = 'options.reference stack is masked out. Change options.mask.t or options.reference stack.';
          error(errormsg);
      end
  end
  
%% Start doing the translation!
fprintf('Starting translation.. ');
% do the translation on all stacks < options.reference_stack and within time mask
% filter reference stack
x = ops.mask.x; y = ops.mask.y; z = ops.mask.z; t = ops.mask.t;
thisfilter = fspecial('gaussian', 3, 1);
refstack = imfilter(instack(x(1):x(2), y(1):y(2), z(1):z(2),ops.reference_stack), thisfilter);
for  stack_idx = ops.reference_stack-1:-1:ops.mask.t(1);
    % apply filter
    thisstack = imfilter(instack(x(1):x(2), y(1):y(2), z(1):z(2),stack_idx), thisfilter);
    [dx, transimg] = register_translation_multires(...
                        single(refstack),...
                        single(thisstack),...
                        t_ops);
    save([temp_path 'dx.' num2str(stack_idx)], 'dx'); % write dx values
    save([temp_path 'timg.' num2str(stack_idx)], 'transimg'); % write transimg values
    t_ops.dx_start = dx;
    fprintf('%d.. ', stack_idx);
end
% do the translation on all stacks > options.reference_stack and within time mask
for  stack_idx = ops.reference_stack+1:1:ops.mask.t(2);
    % apply filter
    thisstack = imfilter(instack(x(1):x(2), y(1):y(2), z(1):z(2),stack_idx), thisfilter);
    
    [dx, transimg] = register_translation_multires(...
                        single(refstack),...
                        single(thisstack),...
                        t_ops);
    save([temp_path 'dx.' num2str(stack_idx)], 'dx'); % write dx values
    save([temp_path 'timg.' num2str(stack_idx)], 'transimg'); % write transimg values
    t_ops.dx_start = dx;
    fprintf('%d.. ', stack_idx);
end
% add values for reference_stack
dx = [0 0 0];
transimg = single(imfilter(instack(x(1):x(2), y(1):y(2),z(1):z(2), ops.reference_stack), thisfilter));
x = ops.mask.x; y = ops.mask.y; z = ops.mask.z; t = ops.mask.t;
save([temp_path 'dx.' num2str(ops.reference_stack)], 'dx');
save([temp_path 'timg.' num2str(ops.reference_stack)], 'transimg');
fprintf('%d.. DONE! %s', ops.reference_stack, datestr(now()));

%% Assemble values in temporal order

  % dx values
  clear dx;
  dx = struct;
  for dx_idx = ops.mask.t(1):ops.mask.t(2)
      thisdx = load([temp_path 'dx.' num2str(dx_idx)], '-mat');
      if ~isfield(dx, 'array')
          dx.array = {thisdx.dx};
      else
          dx.array(end+1) = {thisdx.dx};
      end
      if ~isfield(dx,'indices')
          dx.indices = dx_idx;
      else
          dx.indices(end+1) = dx_idx;
      end
  end
  save([base_path 'dx.mat'], 'dx');
 
  
  % translated images
  fprintf('\nSaving data.. ');
  for img_idx = ops.mask.t(1):ops.mask.t(2)
      tempimg = load([temp_path 'timg.' num2str(img_idx)],'-mat');
      count = fwrite(fid,uint16(tempimg.transimg),'uint16');
      fprintf('%d..',img_idx);
  end
  fprintf('\n%s saved successfully @ %s\n', ops.savefile, datestr(now()));
  % TO DO: automate process of making an altered '.imagine' file for
  % stackmm and other programs to utilize.  FOR NOW, do this manually.

%--------------------------------------------------------------------------