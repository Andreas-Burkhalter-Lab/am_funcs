function ipout = imphysprocess(ipin,ipp)
% IMPHYSPROCESS: apply processing to imphys data
%
% This function can be used to do averaging, background subtraction,
% filtering, and resampling.
%
% Syntax:
%   ipout = imphysprocess(ipin,ipp)
% where
%   ipin is an input IMPHYS structure;
%   ipp is an input IMPHYSP structure, containing the instructions about
%     the nature of the desired processing (see IMPHYSP);
% and
%   ipout is the output IMPHYSP structure containing the result.
%
% See also: IMPHYSP.
  
% Copyright 2005 by Timothy E. Holy
  
  nstacks = length(ipp);
  for i = 1:nstacks
    % For each desired frame, find the input frames which will contribute
    if isfield(ipin,'imfile')
      fileindx = strmatch(ipp(i).oimfile,{ipin.imfile},'exact');
    else
      fileindx = 1:length(ipin);
    end
    stackindx = findainb(ipp(i).ostacknum,[ipin(fileindx).stacknum]);
    stackindx = fileindx(stackindx);
    ipin_tmp = ipin(stackindx);  % These are the frames we need
    % Verify that necessary parameters are identical across contributing
    % frames, and copy those parameters
    ipout(i) = imphyscopy(ipin_tmp,{'dimension'},ipp(i));

    %
    % Frame/stack averaging, subtracting, etc.
    %
    ipout(i).image = ipp(i).ostackweight(1)*single(imphysfetch(ipin_tmp(1)));
    for j = 2:length(ipin_tmp)
      ipout(i).image = ipout(i).image + ...
          ipp(i).ostackweight(j)*single(imphysfetch(ipin_tmp(j)));
    end
    
    %
    % Filtering
    %
    if (isfield(ipp,'filterwidth') & ~isempty(ipp(i).filterwidth))
      if (isfield(ipp,'filtersize') & ~isempty(ipp(i).filtersize))
        filtersize = ipp(i).filtersize;
      else
        filtersize = ipp(i).filterwidth*3;
      end
      smoothfilt = fspecial('gaussian',filtersize,ipp(i).filterwidth);
      ipout(i).image = imfilter(ipout(i).image,smoothfilt);
    end
    
    %
    % Resampling
    %
    if (isfield(ipp,'resampmag') & (ipp(i).resampmag ~= 1))
      mag = ipp(i).resampmag;
      method = 'bilinear';
      if (isfield(ipp,'resampmethod') & ~isempty(ipp(i).resampmethod))
        method = ipp(i).resampmethod;
      end
      ipout(i).image = imresize(ipout(i).image,mag,...
                                method,0);
% $$$       % We have to correct the scales
% $$$       fn = {'x2um','y2um','z2um'};
% $$$       for j = 1:length(fn)
% $$$         if isfield(ipout(i).(fn{j}))
% $$$           ipout(i).(fn{j}) = ipout(i).(fn{j})/mag;
% $$$         end
% $$$       end
% $$$       fn = {'xrange','yrange','zrange'};
% $$$       pos = [2 1 3];
% $$$       for j = 1:length(fn)
% $$$         if isfield(ipout(i).(fn{j}))
% $$$           ipout(i).(fn{j}) = [1 size(ipout(i).image,pos(j))];
% $$$         end
% $$$       end
    end
    
    
    