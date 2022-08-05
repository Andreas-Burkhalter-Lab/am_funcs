function stimuli=parse_valvelabels(labels)
% PARSE_VALVELABELS: extract stimulus details from string
% Syntax:
%   stimuli = parse_valvelabels(labels)
% where
%   labels is a cell array of strings, each describing a particular
%     stimulus
% and
%   stimuli is a cell array of structures, each element consisting of a
%     systematic characterization/parametrization of the stimulus identity.
%
% This function obviously needs to be able to parse your strings properly;
% you are encouraged to use formats that are recognized. In cases where
% this isn't feasible, edit the file "stim_lookup.txt" to add the required
% information to the list.

% Copyright 2007 Zhongsheng Guo
% 2008-11-15: split out as separate function (TEH), small extensions made

   negctrls = {'ringer','vehicle'};  % all lowercase

   if ~iscell(labels)
     labels = {labels};
   end
   
   stimuli={};
   
   for idxLabel=1:length(labels)
      label=labels{idxLabel};
      label=lower(label);
      
      if(isempty(label))
        continue
      end
      
      % remove portion before the last = and = itself
      tpos=find(label=='=');
      if(~isempty(tpos))
         label=strtrim(label(tpos(end)+1:end));
      end
      
      % if there + in the middle, don't know what to do
      tpos=find(label=='+');
      if(~isempty(tpos) && tpos(end)>1)
         stimuli{end+1}=lookup_valve_label(label);
         continue;
      end
      
      % replace w/ with +, and replace w/o with -
      if(strncmp(label, 'w/o', 3))
         label=['-' label(4:end)];
      elseif(strncmp(label, 'w/', 2))
         label=['+' label(3:end)];
      end
            
      % ends with /ddd ?
      tExpression= '\/(\d+)$';
      [s, f, t, match, tokens]=regexp(label, tExpression, 'once');
      if(~isempty(t))
         % found /ddd at the end
         tstim.concentration=1/str2double(tokens{1});
         tstim.conc_unit='relative';
         labels_wo_conc=label(1:t(1)-2);
         stim=lookup_valve_label(labels_wo_conc);

         [stim.concentration]=deal(tstim.concentration);
         [stim.conc_unit]    =deal(tstim.conc_unit);
         stimuli{end+1}=stim;
         continue;
      end
      
      % ends with ddd microM or pM or nM?
      texpr='(\d+)\s*(u|micro|p|n|m)m$';
      [s, f, t, match, tokens]=regexp(label, texpr, 'once');
      if(~isempty(t))
         tstim.conc_unit='M';
         if(isequal(tokens{2}, 'micro') || isequal(tokens{2},'u'))
            tstim.conc_unit='\muM';
         elseif(isequal(tokens{2}, 'p'))
            tstim.conc_unit='pM';
         elseif(isequal(tokens{2}, 'n'))
            tstim.conc_unit='nM';
         elseif(isequal(tokens{2}, 'm'))
            tstim.conc_unit='mM';
         else
            error('coding error');
         end
         tstim.concentration=str2double(tokens{1});

         labels_wo_conc=strtrim(label(1:t(1)-1));
         stim=lookup_valve_label(labels_wo_conc);

         [stim.concentration]=deal(tstim.concentration);
         [stim.conc_unit]    =deal(tstim.conc_unit);
         stimuli{end+1}=stim;
         
         continue;
      end

      % ends w/ 10_dd 10-dd or 10_ddM or 10-ddm or 10_dd m?
      texpr='(10)(?:_|-)(\d+)\s*m?$';
      [s, f, t, match, tokens]=regexp(label, texpr, 'once');
      if(~isempty(t))
         tstim.concentration=10^-str2double(tokens{2});
         tstim.conc_unit='M';
         labels_wo_conc=strtrim(label(1:t(1)-1));
         stim=lookup_valve_label(labels_wo_conc);

         [stim.concentration]=deal(tstim.concentration);
         [stim.conc_unit]    =deal(tstim.conc_unit);
         stimuli{end+1}=stim;
         
         continue;
      end
      
      if any(cellfun(@(s) ~isempty(findstr(s,lower(labels{idxLabel}))),negctrls))
        stim=stimdefault_syn_francesco; stim.category='negative_control';
        stimuli{end+1}=stim;
        continue;
      end
      
      % now lookup anyway
      labels_wo_conc=label;
      stim=lookup_valve_label(labels_wo_conc);
      stimuli{end+1}=stim;
      
   end % for, each valve label
   
function stim=lookup_valve_label(labels_wo_conc)
   % let's see if we can tell the identity
   label=labels_wo_conc;
   
   % fraction number?
   texpr='^(\d+)$';
   [s, f, t, match, tokens]=regexp(label, texpr, 'once');
   if(~isempty(t))
      % yea, frac # 
      stim=frac_num_to_stim(str2double(tokens{1}));
      return; 
   end
   
   % mixture of fractions?
   texpr='^(?:fr)?(\d+)(?:to|-)(\d+)$';
   [from, to, tokenExtents, match, tokens]=regexp(label, texpr, 'once');
   if(~isempty(tokenExtents))
      % yea
      stims={};
      for fra_num=str2double(tokens{1}):str2double(tokens{2})
         stims{end+1}=frac_num_to_stim(fra_num);
      end
      stim=cat(2, stims{:});
      return; 
   end
   
   % now we have to lookup 
   lookupfile=which('stim_lookup.txt');
   
   [status, content]=load_text_file(lookupfile);
   
   expressions=key2value(content, label);
   
   expressions=strtrim(expressions);
   
   if(isempty(expressions))
      error([label ' is not in the lookup table']);
   end
   
   tt=safe_eval([expressions ';'], 'stim');
   
   stim=tt.stim;
   
   stim=fill_and_sort_stim_fields(stim);
   
   
function stim=fill_and_sort_stim_fields(stim)   
   default_fields=struct('category', 'urine_whole', ...
      'identity', '', ...
      'concentration', 1, ...
      'conc_unit', 'relative', ...
      'duration', nan, ...
      'user_tag', '', ...
      'strain', 'Balb/c', ...
      'sex', 'F', ...
      'procedure', {{}} ...
      );
   
   valid_fields=fieldnames(default_fields);
   if(~all(ismember(fieldnames(stim), valid_fields)))
      error('some field are not recognized');
   end
   
   for idxField=1:length(valid_fields)
      fn=valid_fields{idxField};
      if(~isfield(stim, fn))
         [stim.(fn)]=deal(default_fields.(fn));
      end
   end
   
   stim=orderfields(stim, valid_fields);
   
   
   
function stim=frac_num_to_stim(fra_num)
   stim=struct('category', 'urine-derived', ...
      'identity', num2str(fra_num), ...
      'concentration', 1, ...
      'conc_unit', 'relative', ...
      'duration', nan, ...
      'user_tag', '', ...
      'strain', '', ...
      'sex', '', ...
      'procedure', {{}} ...
      );

   stim=fill_and_sort_stim_fields(stim);
   

   
   
   
   
   
   
   
   

