function outValue = keystepper(values,inValue,options)
% KEYSTEPPER: use keypresses to cycle through values
% When looking at a large amount of data, it can be useful to do some
% plotting and then pause for user input before going on to the next
% value.  This function gives the user the ability to navigate among
% values:
%      next value: 'space' or 'right arrow'
%      previous value: 'left arrow'
%      first value: 'up arrow'
%      last value: 'down arrow'
%      arbitrary value: type numbers and press return
%      quit: 'q'
%
% You use this function in the following way:
%    i=1; while ~isnan(i), plot(something(i)); title(num2str(i)); i=keystepper(1:30,i); end
% The 1:30 input specifies the allowed values of i.
%
% Syntax:
%    outValue = keystepper(allowedValues,inValue)
%    outValue = keystepper(allowedValues,inValue,options)
% options may have the following fields:
%   q_only (default false): if true, keystepper returns nan only when the
%     user presses 'q', and stays "pegged" at the end points (without
%     quitting) if the user gets to the beginning or end of the sequence.
%
% See also: GETKEY.
  
% Copyright 2006 by Timothy E. Holy
  
  if (nargin < 3)
    options = struct;
  end
  
  options = default(options,'q_only',false);
  
  currentIndex = find(values == inValue);
  if isempty(currentIndex)
    error('Input value does not fall in acceptable list');
  end
  if (length(currentIndex) > 1)
    error('Non-unique input value');
  end
  callstr = 'set(gcbf,''Userdata'',double(get(gcbf,''Currentcharacter''))) ; uiresume ';
  set(gcf,'keypressfcn',callstr,...
    'userdata','timeout');
  ch = ks_getkey;
  if (ch >= 48 && ch <= 57)
    % User pressed a number, wait for next key press
    digits = char(ch);
    ch = ks_getkey;
    while (ch >= 48 && ch <= 57)
      digits(end+1) = char(ch);
      ch = ks_getkey;
    end
    number = sscanf(digits,'%d');
    cIndex = find(values == number);
    if ~isempty(cIndex)
      currentIndex = cIndex;
    end
  else   
    switch ch
      case 28  % left arrow
        currentIndex = currentIndex-1;
      case {29,32}  % right arrow, space
        currentIndex = currentIndex+1;
      case 30 % up arrow
        currentIndex = 1;
      case 31 % down arrow
        currentIndex = length(values);
      case 112 % p (print fig to temp file)
        tempfilename = ['keystepper_temp_' num2str(currentIndex) '.jpg'];
        print(gcf,'-djpeg100',tempfilename);
        fprintf(['Printed ' tempfilename ' to current directory.\n'])
      case 113 % q
        currentIndex = NaN;
      otherwise
        % currentIndex = NaN;
        % warning('Unrecognized input, exiting keystepper');
        h=show_help([], [], 'keystepper', ...
          [colorize('blue', 'left arrow') ': previous;' 10 ...
          colorize('blue', 'right arrow/space') ': next;' 10 ...
          colorize('blue', 'up arrow') ': first;' 10 ...
          colorize('blue', 'down arrow') ': last;' 10 ...
          colorize('blue', '<type number, press return>') ': go to given value;' 10 ...
          colorize('blue', 'p') ': print current fig to temp .jpg file;' 10 ...
          colorize('red', 'q') ': quit.']);
        uiwait(h);
    end
  end
  if (currentIndex < 1)
    if options.q_only
      currentIndex = 1;
    else
      currentIndex = nan;
    end
  elseif (currentIndex > length(values))
    if options.q_only
      currentIndex = length(values);
    else
      currentIndex = nan;
    end
  end
  if isnan(currentIndex)
    outValue = NaN;
  else
    outValue = values(currentIndex);
  end
end

function ch = ks_getkey
  ch = [];
  while isempty(ch)
    try
      uiwait;
      ch = get(gcf,'UserData');
    catch
      ch = [];
    end
  end
end
