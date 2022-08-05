function h = imreadheader(headerfile)
%  todo: get imrange field

  % This will need to be made more generic to handle .tifs, confocal
  % stacks, etc.
  [s, txt]=load_text_file(headerfile);
  if(s~=0)
    error(['Error opening file: ' headerfile]);
  else
    h.date = key2value(txt,'date and time');
    h.height = str2num(key2value(txt,'image height'));
    h.width = str2num(key2value(txt,'image width'));
    % Find the number of stacks
    strCap=key2value(txt, 'capture params');
    [tok,rem] = strtok(strCap,':;');
    [tok,rem] = strtok(rem,':;');
    h.nstacks = str2num(tok);
    % Frame interval
    [funits,rem] = strtok(rem,':;');
    [fvalue,rem] = strtok(rem,':;');
    fvalue = str2num(fvalue);
    if strmatch(funits,'fps')
      fvalue = 1.0/fvalue;
    elseif ~strmatch(funits,'spf')
      error('Don''t recognize frame interval units');
    end
    h.stacktime = (0:h.nstacks-1)*fvalue;
    % Image file format
    h.machfmt = key2value(txt,'machine format');
    if isempty(h.machfmt)
      h.machfmt = 'n';  % "native"
    end
    bitdepth = key2value(txt,'bit depth');
    if ~isempty(bitdepth)
      bitdepth = str2num(bitdepth);
      if (bitdepth > 8)
        h.prec = 'uint16';
        h.nbytes = 2;
      else
        h.prec = 'uint8';
        h.nbytes = 1;
      end
    else
      h.prec = 'uint16';
      h.nbytes = 2;
    end
    h.filefmt = 'raw';
    h.camera=strtrim(key2value(txt,'camera'));
    strPixelOrder = key2value(txt,'pixel order');
    nativeorder = {'y','x','z'};
    rem = strPixelOrder;
    pixelorder = [];
    while ~isempty(rem)
      [tok,rem] = strtok(rem);
      if(~isempty(tok))
         pixelorder(end+1) = strmatch(tok,nativeorder);
      end % check tok also, in case rem is trailing spaces before call strtok()
    end
    h.pixelorder = pixelorder;
    % Read stimulus info
    strStim = key2value(txt, 'stimulus sequence, time in frames');
    h.stim = eval(['[' strStim ']' ]);
  end
