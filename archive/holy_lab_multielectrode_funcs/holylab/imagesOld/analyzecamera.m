function varargout = analyzecamera(fls,flst,options)
  if (nargin < 2 | isempty(flst))
    flst = cell(1,length(fls));
  end
  if (nargin < 3)
    options = struct;
  end
  fls{:}
  t = input('Enter the exposure times: ');
  if ~isfield(options,'imrange')
    options.imrange = input('Enter the pixel value range (e.g., [0 4095]): ');
  end
  for i = 1:length(fls)
    ipc{i} = imphysfrom2d(fls{i},flst{i});
    ipc{i} = ipc{i}(3:end);  % Toss the first two frames, it can be weird
                             % sometimes
    for j = 1:length(ipc{i})
      ipc{i}(j).inttime = t(i);
      if isfield(options,'imrange')
        ipc{i}(j).imrange = options.imrange;
      end
    end
  end
  ip = cat(2,ipc{:});
  if (length(t) > 1)
    [varargout{1},varargout{2}, varargout{3}] = imnoiseanalysis(ip,options);
  else
    varargout{1} = imnoiseanalysis(ip,options);
  end
  
