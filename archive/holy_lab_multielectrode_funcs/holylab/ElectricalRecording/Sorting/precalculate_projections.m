function precalculate_projections(dirname,op)
% PRECALCULATE_PROJECTIONS: save interactive time in cass
%
% Syntax:
%   precalculate_projections
%   precalculate_projections(dirname)
%   precalculate_proejctions(dirname,options)
%   precalculate_projections(options)
%
% Options: 
%   resumeon = 1 (default) or 0
%       controls whether or not will skip channels that already have a
%       proj.mat file
%   doinorder = 1 (default) or 0
%       reorders the channels to do them in numeric not alphanumeric order
%       (so doesn't eg wait until the very end to do #9)
%
% When executed within a sorting directory (or when passed the dirname of
% a sorting directory), this loops through all the channels and
% pre-computes the data about snippet waveform shape that will be needed
% in CASS, saving it to a file named proj.mat.  This saves interactive
% time.

% Copyright 2007 by Timothy E. Holy

%% input & options

if nargin == 0
    dirname = '.';
    op = struct;
elseif nargin == 1
    if isstruct(dirname)
        op = dirname;
        dirname = '.';
    else
        op = struct;
    end
end
op = default(op,'resumeon',1); % if you need to redo projections, turn this off
op = default(op,'doinorder',1); % slightly slower to start, but means do 1-9 before 10 etc

%% do the work

  dorig = pwd;
  cd(dirname)
  dcur = pwd;
  blocksize = 5000;  % process in blocks of 5000 snippets
  load overview
  sorthead = cass_fix_path(sorthead);
  n_files = length(sorthead);
  channames = dirbyname('chan*');
  if op.doinorder
      channums = NaNs(size(channames));
      for nthChannum = 1:length(channums)
          channums(nthChannum) = str2num(channames{nthChannum}(5:end));
      end
      channums = sort(channums);
      for nthChannum = 1:length(channums)
          channames{nthChannum} = ['chan' num2str(channums(nthChannum))];
      end
  end
  n_chans = length(channames);
  for chanIndex = 1:n_chans
    this_channel = str2num(channames{chanIndex}(5:end));
    this_channel_shIndex = find(this_channel == sorthead(1).channels);
    cd(dcur);
    cd(channames{chanIndex});
    if exist('proj.mat','file') & op.resumeon
      continue
    end
    load autosort_info
    pd = sort_info.projectDirections;
    shc = sortheader_importchan(sorthead,this_channel);
    
    snipProj = {};
    snipmm = {};
    for fileIndex = 1:n_files
      n_snips_per_file(fileIndex) = sorthead(fileIndex).numofsnips(this_channel_shIndex);
    end
%     cum_n_snips_per_file = [0 cumsum(n_snips_per_file)];
%     offset = 1;
    
%     snipIndicesUsed=zeros(2,0); % the used snip indices, top row is file indices

    for fileIndex = 1:n_files
      offset = 1;
      snipProjTmp = [];
      snipmmTmp = [];
      while (offset <= n_snips_per_file(fileIndex))
        rangemax = min(n_snips_per_file(fileIndex),offset+blocksize-1);
        snipIndexTmp = offset:rangemax;
        sniptmp = sortheader_readsnips(shc(fileIndex),snipIndexTmp);
        snipProjTmp{end+1} = pd'*sniptmp;
        snipmmTmp{end+1} = [min(sniptmp); max(sniptmp)];
        offset = offset+blocksize;
      end
      if iscell(snipProjTmp)
        snipProj{fileIndex} = cat(2,snipProjTmp{:}); % cat snipProj from different blocks
        snipmm{fileIndex} = cat(2,snipmmTmp{:});
      else
          snipProj{fileIndex} = snipProjTmp;
          snipmm{fileIndex} = snipmmTmp;
      end
    end
    save proj snipProj snipmm
  end
  cd(dorig)
  
%     while (offset < cum_n_snips_per_file(end))
%       snipIndexTmp = cell(1,n_files);
%       for fileIndex = 1:n_files
%         rangemin = max(1,offset-cum_n_snips_per_file(fileIndex));
%         rangemax = min(n_snips_per_file(fileIndex),...
%           offset+blocksize-1-cum_n_snips_per_file(fileIndex));
%         snipIndexTmp{fileIndex} = rangemin:rangemax;
%        
%         snipIndicesUsed=[snipIndicesUsed ...
%            [fileIndex*ones(1, length(snipIndexTmp{fileIndex})); snipIndexTmp{fileIndex}] ];
%       end
%       sniptmp = sortheader_readsnips(shc,snipIndexTmp);
%       sniptmp = cat(2,sniptmp{:}); % cat sniptmp from different files
%       snipProj{end+1} = pd'*sniptmp;
%       snipmm{end+1} = [min(sniptmp); max(sniptmp)];
%       offset = offset+blocksize;
%     end
%     snipProj = cat(2,snipProj{:}); % cat snipProj from different blocks
%     snipmm = cat(2,snipmm{:});
%     save proj snipProj snipmm
%   end
%   