function base_stacknum = choose_base_stack_interactively(smm,base_stacknum,stacknums,rng)
  % choose_base_stack_interactively: present base stack for visual inspection and approval
  % Syntax:
  %   base_stacknum = choose_base_stack_interactively(smm)
  %   base_stacknum = choose_base_stack_interactively(smm,base_stacknum)
  %   base_stacknum = choose_base_stack_interactively(smm,base_stacknum,stacknums)
  %   base_stacknum = choose_base_stack_interactively(smm,base_stacknum,stacknums,rng)
  % where
  %   smm is the input stackmm object, or 4-D image stack
  %   base_stacknum is the candidate base stack number (default, the middle
  %     of the sequence)
  %   stacknums is a vector of stack numbers that are valid choices for the
  %     base stack, these will be presented to the user if the default is
  %     unacceptable (default is all stacks).
  %   rng is a cell array containing a coordinate range to display
  %     (default, the whole stack)
  % and on output, base_stacknum is the user's choice of base stack (empty
  %   if the user selects 'cancel').
  
  % Copyright 2012 by Timothy E. Holy
  % (snipped out from find_bad_frames, originally written in 2011)
  
  if isa(smm,'stackmm')
    sz = smm.size;
  else
    sz = size(smm);
  end
  sz(end+1:4) = 1;
  if (nargin < 4)
    for i = 1:3
      rng{i} = 1:sz(i);
    end
  end
  if (nargin < 3)
    stacknums = 1:sz(4);
  end
  if (nargin < 2)
    base_stacknum = ceil(sz(4)/2);
  end
  while true
    base_check_stacknum = [0 -1 -2]+base_stacknum;
    base_check_stacknum = base_check_stacknum(base_check_stacknum > 0);
    stk3 = smm(rng{:},base_check_stacknum);
    stk3 = permute(stk3,[1 2 4 3]);
    rgb = double(stk3) / double(max(stk3(:)));
    mplay(rgb)
    %   answer = questdlg('Is the base stack (red channel) OK?');
    answer = input(sprintf('The base is stack %d and is displayed in the red channel. Is it OK? (yes/no/cancel) ',base_stacknum),'s');
    ws = isspace(answer);  % find first non-space character
    answer = answer(find(~ws,1,'first'));
    switch lower(answer)
      case {'','c'}
        base_stacknum = [];
        return
      case 'n'
        [item,ok] = listdlg('PromptString','Choose a different base stack #',...
          'ListString',cellstr(num2str(stacknums)));
        if (ok == 0)
          return
        end
        base_stacknum = stacknums(item);
      case 'y'
        % User is satisfied, we can proceed
        break
    end
  end
end
