function [bioCells comments timeMarkers] = ...
    collect_fake_channel_cass_result(autosortDir, options)

% Read biological cells out of multichannel cass spike sorting output
%
% SYNTAX: [bioCells comments timeMarkers] = ...
%               collect_fake_channel_cass_result(autosortDir, options)
%
% IN: 
%   autosortDir
%   options -   preferredAsiFilenames = a string list containing the saved
%                   names to look first. If you followed the convention that
%                   cass' results are saved as autosort_info.mat, 1.mat,
%                   2.mat, ...(up to 1000.mat), you can safely ignore this
%                   option.
% 
%               easy_read_comments = 0 (default) or 1; takes the comment
%                   timeMarker information out of the cell format for quick
%                   and easy hand browsing
%               salvage_from_trash = 0 (default) or 1; if 1, a cell
%                   containing some templates that have been merged with
%                   trash and some that haven't will be preserved, with a
%                   "T" for trash added to the comments field of the one(s)
%                   that was/were trashed (not necessarily recommended -
%                   use at your own risk!)
%
% OUT:  
%   bioCells    a cell array with one entry per biological cell identified;
%               each entry contains a list of the templates that are a part
%               of that cell
%   comments    a cell array indexed like bioCells with each entry
%               its own cell array of all comments fields associated
%               with that biological cell; OR, if the easy_read_comments
%               option set, will just give a separated string per cell
%   timeMarkers like comments but for timeMarker info
%
% HISTORY:
%   ?           (ZG)    wrote it
%   2007-03-15  (RCH)   added in comment and timemarker compilation feature
%   ...         (ZG)    please see the svn log.

if (nargin < 2)
  options = struct;
end

if(nargin==0)
   autosortDir='sort_fake'; % the dir that holds autosort result of fake channels
end

% It would seem that what follows could be more elegantly & flexibly done using
% prepare_chanfile---discuss
if(~isfield(options, 'preferredAsiFilenames'))
   options.preferredAsiFilenames={}; % Asi: AutoSort Info. A name appears earlier in the list has higher priority.
   for i=1000:-1:1
      options.preferredAsiFilenames{end+1}=[num2str(i) '.mat'];
   end
end

if(~iscell(options.preferredAsiFilenames))
   options.preferredAsiFilenames={options.preferredAsiFilenames};
end
if(isempty(strmatch('autosort_info.mat', options.preferredAsiFilenames, 'exact')))
   options.preferredAsiFilenames{end+1}='autosort_info.mat';
end

options = default(options,'salvage_from_trash',false);
options = default(options,'easy_read_comments',false);

overviewFile=fullfile(autosortDir, 'overview.mat');
overview=load(overviewFile, '-mat');
sorthead=overview.sorthead(1); % the first ssnp is enough b/c they share the same template file
sorthead = cass_fix_path(sorthead,options);
ssnpFile=fullfile(sorthead.fh.abspathstr, sorthead.fh.filename);
fakeChannelInfoFile=replace_extension(ssnpFile, '.fake_info');
tt=load(fakeChannelInfoFile, '-mat');
fakeInfo=tt.fakeInfo;
templateFile=fullfile(sorthead.fh.abspathstr,tt.templateFile);

templateVars=load(templateFile, '-mat');
templates=templateVars.templates;
nTemplates=length(templates);

nChannelsToProcess=nTemplates;
% nChannelsToProcess=1; % for debug purpose

% mergeMatrix=zeros(nTemplates, nTemplates);
mergeMatrix=eye(nTemplates);
comments_holder = cell(1,nTemplates);
timeMarkers_holder = cell(1,nTemplates);
isTrash=zeros(1, nTemplates);
trashedStars=cell(1, nTemplates); % a cell array of vectors, which holds the reason why a template was trashed.
for templateID=1:nTemplates
   trashedStars{templateID}=[];   
end

for fakeChannel=1:nChannelsToProcess
   templateID=fakeChannel;
   
   autosortInfoFileParentDir=fullfile(autosortDir, ['chan' num2str(fakeChannel)]);
   autosortInfoFile=find_first_pattern(options.preferredAsiFilenames, ...
      autosortInfoFileParentDir);
   nFiles=length(autosortInfoFile);
   if(nFiles==0)
      error('can not find auto sort info file');
   elseif(nFiles>1)
      error('found multiple candidates for auto sort info file');
   end
   autosortInfoFile=autosortInfoFile{1};
   
   tt=load(autosortInfoFile, '-mat');
   landmarkClust=tt.sort_info.landmarkClust;
   if isfield(tt.sort_info,'comment')
   comments_holder{fakeChannel} = tt.sort_info.comment;
   end
   timeMarkers_holder{fakeChannel} = tt.sort_info.timeMarker;
   
   isMergedWithCurChannel=landmarkClust==landmarkClust(1);
   
   mergedTemplates=fakeInfo(fakeChannel).closestTemplateIDs(isMergedWithCurChannel);
   
   if(mergedTemplates(1)~=templateID)
      error('coding error');
   end
   
   if(landmarkClust(1)==0)
      % merged w/ trash
      isTrash(mergedTemplates)=1;
      for tt=make_vector(mergedTemplates, 'row')
         trashedStars{tt}(end+1)=templateID;
      end
   else
      mergeMatrix(mergedTemplates, templateID)=1;
   end
end % for, each fake channel

% tt=mergeMatrix | mergeMatrix';
% bioCells=connected_components(tt);
bioCells=connected_components(mergeMatrix);

indicesBioCellsToRemove=[];
for idxCell=1:length(bioCells)
    mergedTemplates=bioCells{idxCell};
    if(any(isTrash(mergedTemplates)))
        if ~options.salvage_from_trash 
            indicesBioCellsToRemove(end+1)=idxCell;
            if(length(mergedTemplates)>1)
                trashedTemplates=mergedTemplates(logical(isTrash(mergedTemplates)));
                msg=['templates ' num2str(mergedTemplates) ' are merged as desired cell, \nbut' ...
                    ' some of them was/were also merged with trash;\n' ...
                    'the trashed template(s) is/are: ' num2str(trashedTemplates) '\n'];
                for tt=make_vector(trashedTemplates, 'row')
                   msg=[msg '   template ' num2str(tt) ' is trashed due to ' ...
                      num2str(trashedStars{tt}) '\n'];
                end
                warning(sprintf(msg));
            end
        elseif all(isTrash(mergedTemplates))
            indicesBioCellsToRemove(end+1)=idxCell;
        else
            iTrashed = find(isTrash(mergedTemplates));
            for nthTrashed = 1:length(iTrashed)
                prevComment = comments_holder{mergedTemplates(nthTrashed)};
                if iscell(prevComment)
                    prevComment = prevComment{1};
                end
                comments_holder{mergedTemplates(nthTrashed)} = [prevComment ' T'];
            end
        end
    end
end

bioCells(indicesBioCellsToRemove)=[];
comments = cell(1,length(bioCells));
timeMarkers = cell(1,length(bioCells));
for nthCell = 1:length(bioCells)
    if options.easy_read_comments
        comm = '';
        for nthComm = 1:length(bioCells{nthCell})
            newcomm = comments_holder{bioCells{nthCell}(nthComm)};
            if iscell(newcomm)
                newcomm = newcomm{1};
            end
            comm = [comm newcomm ' /'];
        end
        comments{nthCell} = comm(1:(end-1));
    else
    comments{nthCell} = {comments_holder{bioCells{nthCell}}};
    end
    timeMarkers{nthCell} = {timeMarkers_holder{bioCells{nthCell}}};
end
