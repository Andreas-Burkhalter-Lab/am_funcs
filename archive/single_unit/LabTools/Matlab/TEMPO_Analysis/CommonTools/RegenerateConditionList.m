%   RegenerateConditionList - Sifts through stimulus parameters to regenerate the condition list.
%       Beginning with condition 0, assigns conditions numbers to each unique condition
%       BJP 8/9/01
%       
%       conditions is a an array containing the parameter value for each trial and also includes a row for the condition #
%       UniqueConds is an array containing the corresponding parameter values for each unique condition
%       Param is a cell containing the strings describing each parameter varied
%       num_conditions is a integer indicating the number of unique conditions
function [conditions, UniqueConds, param, num_conditions] = RegenerateConditionList(data, BegTrial, EndTrial, PATH, FILE, Protocol)

TEMPO_Defs;
ProtocolDefs;

num_modality = 0;

% get unique conditions  
if (Protocol < 200 | Protocol > 399)  % >200 is not just BJP anymore.  JWN 021105
    for modality = 1:NUM_DOTS_PARAMS;
   	   % get values for parameter 
	    mod_values = data.dots_params(modality,:,PATCH1);
        mod_values = mod_values(~isnan(mod_values))';
%         ['~isempty(mod_values) equals ' num2str(~isempty(mod_values))]
%         ['size(munique(mod_values,1))' num2str(size(munique(mod_values)))]
%         ['modality equals ' num2str(modality)]
%         ['DOTS_BIN_CORR_SEED' num2str(DOTS_BIN_CORR_SEED)]
        if ((~isempty(mod_values)))
            [a,b] = munique(mod_values);
            [c,d] = size(a);
           ['c' num2str(c)]
           modality ~= DOTS_BIN_CORR_SEED
   	    if ((c ~= 1) & (modality ~= DOTS_BIN_CORR_SEED) )
      	    num_modality = num_modality + 1;
            conditions(num_modality,:) = mod_values'; 
            param{num_modality} = [dots_keywords{modality} ' '];
            param{num_modality} = param{num_modality}(6:end);       % remove 'DOTS_'
            param{num_modality}(find (param{num_modality} =='_')) = ' '; %remove underscore
        end 
        end
    end 
    UniqueConds = munique(conditions');
    UniqueConds = sortrows(UniqueConds,1);
    %start with condition 0
    UniqueConds (:,size(UniqueConds,2) + 1) = [1:size(UniqueConds,1)]' - 1;
    
    num_conditions = size(UniqueConds,1);

    for cond = 1:(num_conditions)
        select_trials = [];
        mod = 1;
        while mod <= num_modality
            %accumulate trials corresponding with nested params
            select_trials = [select_trials find(conditions(mod,:) == UniqueConds(cond, mod))  ];
            mod = mod + 1;    
        end
        conditions (num_modality + 1,select_trials) = cond - 1;
    end
    
       UniqueConds = sortrows(UniqueConds,num_modality + 1);

else
    switch Protocol
    
    case FIX_1_23_45
        param{1} = ' '; %obj_keywords{OBJ_STATUS};
        param{2} = ' '; %obj_keywords{OBJ_STATUS};
        param{3} = ' '; %obj_keywords{OBJ_STATUS};
        param{4} = ' '; %obj_keywords{OBJ_STATUS};
        param{5} = ' '; %obj_keywords{OBJ_STATUS};
        param{6} = 'Background Color '; 
       
        nest_var(1) = OBJ_STATUS;
        nest_var(2) = OBJ_STATUS;
        nest_var(3) = OBJ_STATUS;
        nest_var(4) = OBJ_STATUS;
        nest_var(5) = OBJ_STATUS;
        nest_var(6) = BKGND_BACK_COLOR;
        nest_obj(1) = 1;
        nest_obj(2) = 2;
        nest_obj(3) = 3;
        nest_obj(4) = 4;
        nest_obj(5) = 5;
        nest_obj(6) = 0;    %background parameters varied
        
    case FIX_1_23
        param{1} = ' '; %obj_keywords{OBJ_STATUS};
        param{2} = ' '; %obj_keywords{OBJ_STATUS};
        param{3} = ' '; %obj_keywords{OBJ_STATUS};
        param{4} = 'Background Color '; 
        nest_var(1) = OBJ_STATUS;
        nest_var(2) = OBJ_STATUS;
        nest_var(3) = OBJ_STATUS;
        nest_var(4) = BKGND_BACK_COLOR;
        nest_obj(1) = 1;
        nest_obj(2) = 2;
        nest_obj(3) = 3;
        nest_obj(4) = 0;
        
    case BIND_DIR_TUNING
        param{1} = 'Orientation '; %obj_keywords{OBJ_ROTATION};
        param{2} = 'Direction ' ; %obj_keywords{OBJ_DIR};
        nest_var(1) = OBJ_ROTATION;
        nest_var(2) = OBJ_DIR;
        nest_obj = 1:2;
    case BIND_SPATFREQ_TUNING
        param{1} = 'Spat Freq '; %bj_keywords{OBJ_AMPL};
        nest_var(1) = OBJ_AMPL;
        nest_obj = 1:2;
    case BIND_TEMPFREQ_TUNING
        param{1} = 'Temp Freq '; %obj_keywords{OBJ_FREQ};
        nest_var(1) = OBJ_FREQ;
        nest_obj = 1:2;
    case BIND_HDISP_TUNING
        param{1} = 'H. Disp '; %obj_keywords{OBJ_FREQ};
        nest_var(1) = OBJ_HDISP;
        nest_obj = 1:2;
    case BIND_BAR_DIR_TUNING
        param{1} = 'Site '; %obj_keywords{OBJ_ROTATION};
        param{2} = 'Direction ' ; %obj_keywords{OBJ_DIR};
        nest_var(1) = OBJ_STATUS;
        nest_var(2) = OBJ_ROTATION;
        nest_obj(1) = 1;
        nest_obj(2) = 2;
    case BIND_RF_MAPPING
        param{1} = 'X-CTR '; %obj_keywords{OBJ_ROTATION};
        param{2} = 'Y-CTR ' ; %obj_keywords{OBJ_DIR};
        nest_var(1) = OBJ_TRANS_X;
        nest_var(2) = OBJ_TRANS_Y;
        nest_obj(1) = 1;
        nest_obj(2) = 1;
        param{3} = 'X-CTR '; %obj_keywords{OBJ_ROTATION};
        param{4} = 'Y-CTR ' ; %obj_keywords{OBJ_DIR};
        nest_var(3) = OBJ_TRANS_X;
        nest_var(4) = OBJ_TRANS_Y;
        nest_obj(3) = 2;
        nest_obj(4) = 2;
        
        
        
    case FIX_VARY_HISTORY
        param{1} = ' ';
        param{2} = ' ';
        param{3} = ' ';
        param{4} = ' ';
        param{5} = ' ';
        param{6} = ' ';
        nest_var(1) = OBJ_STATUS;
        nest_var(2) = OBJ_STATUS;
        nest_var(3) = OBJ_STATUS;
        nest_var(4) = OBJ_STATUS;
        nest_var(5) = OBJ_STATUS;        
        nest_var(6) = BKGND_BACK_COLOR;
        nest_obj(1) = 1;
        nest_obj(2) = 2;
        nest_obj(3) = 3;
        nest_obj(4) = 4;
        nest_obj(5) = 5;
        nest_obj(6) = 0;
    case FIX_VARY_BACKCOLOR
        param{1} = ' ';
        param{2} = ' ';
        param{3} = ' ';
        param{4} = ' ';
        nest_var(1) = OBJ_STATUS;
        nest_var(2) = OBJ_STATUS;
        nest_var(3) = OBJ_STATUS;
        nest_var(4) = BKGND_BACK_COLOR;
        nest_obj(1) = 1;
        nest_obj(2) = 2;
        nest_obj(3) = 3;
        nest_obj(4) = 0;   
        
    case FIX_TRANS_QUAD
        param{1} = 'Dir: ';
        param{2} = 'Backcolor: ';
        param{3} = 'Jitter: ';
        nest_var(1) = OBJ_DIR;
        nest_var(2) = BKGND_BACK_COLOR;
        nest_var(3) = OBJ_JITTER;
        nest_obj(1) = 1;
        nest_obj(2) = 0;
        nest_obj(3) = 1;
                
    case BIND_COHER_DISC
        param{1} = 'Direction ' ;
        param{2} = 'Phase Jitter';
        nest_var(1) = OBJ_STATUS;
        nest_var(2) = OBJ_STATUS;
        nest_var(3) = OBJ_STATUS;
        nest_var(4) = OBJ_STATUS;
        nest_var(5) = OBJ_DIR;
        nest_var(6) = OBJ_JITTER;
        nest_obj(1) = 1;
        nest_obj(2) = 2; 
        nest_obj(3) = 3;
        nest_obj(4) = 4; 
        nest_obj(5) = -1;   %-1 means only grab value for object currently on
        nest_obj(6) = -1;   %-1 means only grab value for object currently on
    end    
    
       for modality = 1:length(nest_var)
 		  	% get values for parameter 
            if (nest_obj(modality) == -1)
                %pull values for protocols in which object status is varied, with only single obj on.
                conditions (modality,:) = data.obj_params( OBJ_STATUS ,:,1) .*data.obj_params( nest_var(modality) ,:,1) + data.obj_params( OBJ_STATUS ,:,2) .*data.obj_params( nest_var(modality) ,:,2) + data.obj_params( OBJ_STATUS ,:,3) .*data.obj_params( nest_var(modality) ,:,3) + data.obj_params( OBJ_STATUS ,:,4) .*data.obj_params( nest_var(modality) ,:,4) + data.obj_params( OBJ_STATUS ,:,5) .*data.obj_params( nest_var(modality) ,:,5) + data.obj_params( OBJ_STATUS ,:,6) .*data.obj_params( nest_var(modality) ,:,6);
            elseif (nest_obj(modality) == 0)
                %pull value for nested bkgnd parameters
                conditions (modality,:) = data.bkgnd_params( nest_var(modality) ,:);
            else    
                %pull value for nested obj parameters, multiple objects on
                conditions (modality,:) = data.obj_params( nest_var(modality) ,:,nest_obj(modality) );
            end    
       end
   if ( size(data.misc_params,1) >= CONDITION )
       conditions (end + 1,:) = data.misc_params(CONDITION, :);
   else
        conditions (end + 1,:) = 0;
        if (Protocol == FIX_1_23_45)
            for trial = 1:size(conditions,2)
                if conditions(1:end-2, trial) == [0 0 0 0 0]'
                   conditions(end,trial) = 1;
                end
                if conditions(1:end-2, trial) == [1 0 0 0 0]'
                   conditions(end,trial) = 2;
                end
                if conditions(1:end-2, trial) == [0 1 0 0 0]'
                   conditions(end,trial) = 3;
                end
                if conditions(1:end-2, trial) == [0 0 1 0 0]'
                   conditions(end,trial) = 4;
                end
                if conditions(1:end-2, trial) == [0 1 1 0 0]'
                   conditions(end,trial) = 5;
                end
                if conditions(1:end-2, trial) == [0 0 0 1 0]'
                   conditions(end,trial) = 6;
                end
                if conditions(1:end-2, trial) == [0 0 0 0 1]'
                   conditions(end,trial) = 7;
                end
                if conditions(1:end-2, trial) == [0 0 0 1 1]'
                   conditions(end,trial) = 8;
                end
            end
        end
        if (Protocol == FIX_1_23)
            for trial = 1:size(conditions,2)
                if conditions(1:end-2, trial) == [0 0 0]'
                   conditions(end,trial) = 1;
                end
                if conditions(1:end-2, trial) == [1 0 0]'
                   conditions(end,trial) = 2;
                end
                if conditions(1:end-2, trial) == [0 1 0]'
                   conditions(end,trial) = 3;
                end
                if conditions(1:end-2, trial) == [0 0 1]'
                   conditions(end,trial) = 4;
                end
                if conditions(1:end-2, trial) == [0 1 1]'
                   conditions(end,trial) = 5;
                end
            end
        end
    
    end
   
   if  (isempty(find (conditions(1:end-1,conditions(end,:) == 0) == data.one_time_params(NULL_VALUE))  ) & ( (Protocol == BIND_DIR_TUNING) | (Protocol == BIND_SPATFREQ_TUNING) | (Protocol == BIND_TEMPFREQ_TUNING) | (Protocol == BIND_HDISP_TUNING) ) )
%        UniqueConds(UniqueConds(:,end) == 0,1:end-1) = data.one_time_params(NULL_VALUE);
        conditions(1:end-1,conditions(end,:) == 0) = data.one_time_params(NULL_VALUE);
    end
   num_modality = length(nest_var);
   UniqueConds = munique(conditions');
   UniqueConds = sortrows(UniqueConds,num_modality + 1);
   % some runs did not insert null values into params during control trials
  
   % once have derived which conditions vary, get list of unique stimulus conditions, sort, and get total number
   num_conditions = size(UniqueConds,1);

end


