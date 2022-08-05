function [x,memm] = subsref(memm,s)
% MERECMM/SUBSREF: memory-mapped MEREC files
% Access methods:
%    x = memm(13,:)
%    h = memm.header
%    sz = memm.size
%    type = memm.type (e.g., 'uint16')
  
% Copyright 2006 by Timothy E. Holy
  
  switch s.type
   case '()'
    % We're asking for data
    indx = s.subs;
    switch length(indx)
     case 1
       % Addressing using indexing rather than channel number
       error('Not yet implemented');
       %x = memm.mmla(indx{1});
       %error('Must address by channel, scan#');
      case 2
        % Addressing using channel #, scan #
        % Process the channel information
        if (ischar(indx{1}) && strmatch(indx{1},':'))
          indx{1} = 1:length(memm.header.channels);
        else
          indx{1} = memm.chan2ind(indx{1}+1);  % +1 to account for chan0
          if any(isnan(indx{1}))
            error('Channel not recorded');
          end
        end
        % Process the scan information
        contiguous = memm.contiguous;
        if ischar(indx{2})
          if strmatch(indx{2},':')
            indx{2} = [1 memm.size(2)];
            contiguous = true;
          else
            error('Invalid second character argument');
          end
        end
        if contiguous && length(indx{2}) > 2
          error(['In contiguous reads, second argument must be a 1- or' ...
            ' 2-vector']);
        end
        indx2 = double(indx{2})-1;  % go to zero-offset values
        if contiguous
          nscans_to_read = indx2(end)-indx2(1)+1; % works for scalars & vectors
        else
          nscans_to_read = length(indx2);
        end
        min_indx2 = min(indx2);
        max_indx2 = max(indx2);
        if (min_indx2 < 0 || max_indx2 > memm.size(2))
            error('Invalid range of scan numbers requested');
        end
        blockstart = min_indx2:memm.blocksz_in_scans:max_indx2;
        n_blocks = length(blockstart);
        if ~contiguous && n_blocks > 1
          % Break indx2 into blocks
          [indx2s,sortIndex] = sort([blockstart(:); indx2(:)]);
          indxbreak = [find(sortIndex <= length(blockstart)); ...
            length(indx2s)+1];
          indx2_block = cell(1,n_blocks);
          xindx_block = cell(1,n_blocks);
          for blockIndex = 1:n_blocks
            rng = indxbreak(blockIndex)+1:indxbreak(blockIndex+1)-1;
            indx2_block{blockIndex} = indx2s(rng);
            xindx_block{blockIndex} = sortIndex(rng)-n_blocks;
          end
        end
        if n_blocks > 1
          % Preallocate
          x = zeros(length(indx{1}),nscans_to_read,memm.type);
        end
        for blockIndex = 1:n_blocks
          cstart_in_scans = blockstart(blockIndex);
          cend_in_scans = min(cstart_in_scans+memm.blocksz_in_scans-1,...
            max_indx2);
          % Determine whether this range of scan numbers lies within the
          % mappable region; if not, we have to move the offset
          coffset_in_scans = (memm.mm.offset-memm.offset)/(sizeof(memm.type)*memm.size(1));
          if (cstart_in_scans < coffset_in_scans || ...
              cend_in_scans >= coffset_in_scans + memm.blocksz_in_scans)
            coffset_in_scans = cstart_in_scans;
            memm.mm.offset = coffset_in_scans*memm.size(1)*sizeof(memm.type)+memm.offset;
          end
          % Check to see if the mapped region extends beyond the end of the
          % file; if so, we need to shrink it
          bufsize_in_scans = memm.mm.Format{2}(2);
          if (coffset_in_scans + bufsize_in_scans > memm.size(2))
            bufsize_in_scans = memm.size(2) - coffset_in_scans;
            memm.mm.Format{2}(2) = bufsize_in_scans;
          end
          % A previous read near the end of the file might have truncated
          % the mapped buffer size; check to see if this needs restoring
          nscans_to_read = cend_in_scans - cstart_in_scans + 1;
          if (nscans_to_read > bufsize_in_scans)
            memm.mm.Format{2}(2) = nscans_to_read;
          end
          if contiguous
            orng = [1, nscans_to_read]+cstart_in_scans-coffset_in_scans;
            orng = orng(1):orng(2);
          else
            if (n_blocks > 1)
              orng = indx2_block{blockIndex}-(coffset_in_scans-1);
            else
              orng = indx2-(coffset_in_scans-1);
            end
          end
          % Read the data
          if (n_blocks > 1)
            % We have multiple blocks, so need to index on x, too
            if contiguous
              xrng = [cstart_in_scans cend_in_scans] - indx2(1) + 1;
              xrng = xrng(1):xrng(2);
            else
              xrng = xindx_block{blockIndex};
            end
            x(:,xrng) = memm.mm.Data.x(indx{1},orng);
          else
            % Since there is only one block, we can read in one go
            x = memm.mm.Data.x(indx{1},orng);
          end
        end
      otherwise
        error('Index merecmm with one or two arguments only')
    end
    % Convert to volts if necessary
    if memm.tovolts
      x = (single(x)-memm.header.minSample) * memm.d2v_slope + ...
        memm.d2v_offset;
    end
    case '.'
      switch s.subs
        case 'size'
          x = memm.size;
        case 'type'
          x = memm.type;
        case 'filename'
          x = memm.filename;
        case 'pathstr'
          x = memm.pathstr;
        case 'basestr'
          x = memm.basestr;
        case 'extstr'
          x = memm.extstr;
        case 'blocksize'
          x = memm.blocksz_in_scans * memm.size(1);
        case 'blocksz_in_scans'
          x = memm.blocksz_in_scans;
        case 'tovolts'
          x = memm.tovolts;
        case 'header'
          x = memm.header;
        case 'channels'
          x = memm.header.channels;
        case 'scanrate'
          x = memm.header.scanrate;
        case 'contiguous'
          x = memm.contiguous;
        case 'nscans'
          x = memm.header.nscans;
        otherwise
          error(['Field ' s.subs{1} ' not recognized']);
      end
  end
