function [data comments timeMarkers] = ...
                        load_fake_channel_celltimes(partdirs,options)

% LOAD_FAKE_CHANNEL_CELLTIMES: import multichannel sorting into ephys
%
% Syntax:
%   [data comments timeMarkers] = load_fake_channel_celltimes(partdirs,options)
% where
%   merecfiles is a cell array of the .merec files for this experiment;
%   partdirs is a cell array of directory names for separately-sorted
%     groups of channels;
%   options: see code
% and
%   data is an ephys structure
%   comments is a cell array (one for each biocell) of cell arrays (one
%     entry for each template in the biocell) of the comment fields from
%     cass
%   timeMarkers is structured like comments, but with timeMarker informatio
%
% HISTORY:
%   ?           (TH)    wrote it
%   2007-04-14  (RCH)   added ability to recognize merec basenames modified
%                       with _part1_fake as well as _part2_fake; added
%                       ability to find .fit files that have been
%                       sequestered in a "fitted" directory as per
%                       francesco's arrangement; added salvage_from_trash
%                       option; added ability to output timemarkers and comments
%   2007-04024  (RCH)   resized the celltimes field to conform to the
%                       output that's been standard for cass until now
%
% See also: vnosorted2ephys, vno2dcell
  
% Copyright 2007 by Timothy E. Holy

  
  if (nargin < 2)
    options = struct;
  end
  options = default(options,'tmin',0.002);  % 2ms double-trigger protection
  options = default(options,'salvage_from_trash',0); 
    % turning this on means that clusters that were in one place saved as
    % good an in another place merged with trash won't result in everything
    % they're a part of being trashed; instead, the trashed status will be
    % indicated in the comments field

  overviewFile=fullfile(partdirs{1}, 'overview.mat');
  overview=load(overviewFile, '-mat');
  sorthead = cass_fix_path(overview.sorthead);
  merecfiles=getMerecfiles(sorthead);
  
  % get scanrate
  for fileIndex=1:length(merecfiles)
     header=readheader(merecfiles{fileIndex});
     scanrates(fileIndex)=header.scanrate;
  end
  
  data = ephysfromai(merecfiles);
  data = rmfield(data,{'snipfile','cellnums'});
%   if isfield(data,'envelopefile')
%     data = rmfield(data,'envelopefile');
%   end
  [data.sort_template] = deal(true);
  
  for pdIndex = 1:length(partdirs)
    [mccells comments{pdIndex} timeMarkers{pdIndex}] = ...
        collect_fake_channel_cass_result(partdirs{pdIndex},options);
    n_cells = length(mccells);
    overviewFile=fullfile(partdirs{pdIndex}, 'overview.mat');
    overview=load(overviewFile, '-mat');
    sorthead = cass_fix_path(overview.sorthead);
    n_files = length(sorthead);
    
    newMerecfiles=getMerecfiles(sorthead);
    if(~isequal(merecfiles, newMerecfiles))
       error('Different parts used different merec file set');
    end
    
    % preload all sniptimes
    fitResultFiles=getFitResultFiles(sorthead);
    templateFile=getTemplateFile(fitResultFiles{1});
    tt=load(templateFile, '-mat');
    channels=tt.channels;
    fineClusters=tt.fineClusters;
    nTotalTemplates=length(tt.templates);
    sniptimes=cell(n_files, nTotalTemplates);
    for fileIndex=1:n_files
       newTemplateFile=getTemplateFile(fitResultFiles{fileIndex});
       if(~isequal(newTemplateFile, templateFile))
          error('More than one template files are used for the .fit file set'); 
       end
       fprintf('Loading fitting result from %s ...', fitResultFiles{fileIndex});
       if ~exist(fitResultFiles{fileIndex},'file')
           fitResultFiles{fileIndex} = ['fitted' filesep fitResultFiles{fileIndex}];
       end
       tt=load(fitResultFiles{fileIndex}, 'fitting', '-mat');
       fprintf(' done.\n');
       fitting=tt.fitting;
       for templateID=1:nTotalTemplates
          sniptimes{fileIndex, templateID}=[fitting{templateID}.shiftedTime];
       end % for, each template
    end % for each .fit file
    
    % now load primary channels
    primaryChannels=cell(1,nTotalTemplates);
    for templateID=1:nTotalTemplates
       primaryChannels{templateID}=channels(fineClusters(templateID).channelsUsedToAlign);
    end % for, each templateID
    
    for cellIndex = 1:n_cells
      n_templates = length(mccells{cellIndex});
      t_agg = cell(1,n_files);
      curCellPrimChans=cell(1,n_templates);
      for templateIndex = 1:n_templates
        % For each template that contributes to the current cell, load
        % the spike times
        %shc = sortheader_importchan(sorthead, ...
        %  mccells{cellIndex}(templateIndex));
        templateID=mccells{cellIndex}(templateIndex);
        
        % fill in primary channel info for cur cell
        curCellPrimChans{templateIndex}=primaryChannels{templateID};
        
        % Aggregate the times with the times from the previous template
        % that is contributing to this cell
        for fileIndex = 1:n_files
           cur_sniptimes = sniptimes{fileIndex, templateID};
           t_agg{fileIndex} = [t_agg{fileIndex} ...
              cur_sniptimes];
              %double(shc(fileIndex).sniptimes)];
        end
      end
      for fileIndex = 1:n_files
        % We have to sort the times because they came from different
        % templates
        t_agg{fileIndex} = sort(t_agg{fileIndex});
        % discard double triggers
        dt = diff(t_agg{fileIndex});
        killIndex = find(dt < options.tmin*scanrates(fileIndex))+1;
        t_agg{fileIndex}(killIndex) = [];
      end
      % Give the cell a name
      celltag = [partdirs{pdIndex} sprintf(';%d',mccells{cellIndex})];
      if isfield(data,'cellnums')
        cellnum = length(data(1).cellnums)+1;
      else
        data(1).cellnums = [];
        cellnum = 1;
      end
      % Store the cell data
      for fileIndex = 1:n_files
        data(fileIndex).cellnums(end+1) = cellnum;
        if ~isfield(data,'celltags')
          data(fileIndex).celltags = {celltag};
        else
          data(fileIndex).celltags{end+1} = celltag;
        end
        
%         if(~isfield(data, 'primaryChannels'))
%            data(fileIndex).primaryChannels{1}=curCellPrimChans; 
%         else
%            data(fileIndex).primaryChannels{end+1}=curCellPrimChans;
%         end
        if(~isfield(data, 'cellchandef'))
           data(fileIndex).cellchandef{1}=cat(2,curCellPrimChans{:}); 
        else
           data(fileIndex).cellchandef{end+1}=cat(2,curCellPrimChans{:});
        end
        
        % fix cellchandef
        %data(fileIndex).cellchandef{end+1} =
        if ~isfield(data,'celltimes')
          data(fileIndex).celltimes{1} = t_agg{fileIndex}';
        else
          data(fileIndex).celltimes{end+1} = t_agg{fileIndex}'; % t_agg{i} is the cur cell's celltimes for files{i}
        end
      end % for, each file
      
    end % for, each cell
    
  end % for, each part directory
  
  sortedMerecs=sort_ai_by_time(merecfiles);
  indices=findainb(sortedMerecs, merecfiles);
  data=data(indices);
  comments = cat(2,comments{:});
  timeMarkers = cat(2,timeMarkers{:});
  

function merecfiles=getMerecfiles(sorthead)
  merecfiles=cell(1, length(sorthead));
  for idx=1:length(sorthead)
     fh=sorthead(idx).fh;
     if(string_endswith(fh.filename, '_part1_fake.ssnp')) | ...
             (string_endswith(fh.filename, '_part2_fake.ssnp')) 
        filename=[fh.filename(1:end-length('_part2_fake.ssnp'))  '.merec'];
     elseif(string_endswith(fh.filename, '_fake.ssnp')) 
        filename=[fh.filename(1:end-length('_fake.ssnp'))  '.merec'];
     else
        error('can not guess merec file name');
     end
     merec=find_first_pattern({fullfile(fh.abspathstr, filename), ...
        fullfile(fh.abspathstr, '/../', filename)});
     if(isempty(merec))
        error('can not find a merec file');   
     end
     merecfiles{idx}=canonicalize_file_name(merec{1});
  end
  
function fitResultFiles=getFitResultFiles(sorthead)
  fitResultFiles=cell(1, length(sorthead));
  for idx=1:length(sorthead)
     fh=sorthead(idx).fh;
     if(string_endswith(fh.filename, '_fake.ssnp'))
        filename=[fh.filename(1:end-length('_fake.ssnp')) '.fit'];
     else
        error('can not guess .fit file name');
     end
     fitResultFiles{idx}=fullfile(fh.abspathstr, filename);
  end

function templateFile=getTemplateFile(fitResultFile)
   if ~exist(fitResultFile,'file')
       fitResultFile = ['fitted' filesep fitResultFile];
   end
   tt=load( fitResultFile, 'templateFile', '-mat');
   templateFile=replace_filename(fitResultFile, tt.templateFile);
   