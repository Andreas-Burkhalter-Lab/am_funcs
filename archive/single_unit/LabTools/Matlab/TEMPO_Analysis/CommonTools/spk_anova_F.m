function [P_value, MSgroup, MSerror] = spk_anova(data_cnts, trials_group, category_group)

% SPK_ANOVA Based on the values of a parameter (category group, e.g. h disp),
%				takes trials with particular parameters (trials_group) and spike 
%				counts (data_cnts), carries out a one-way ANOVA, suppresses
%				figure output of the results and returns the P value (P_value).  
%			 

temp_matrix = trials_group;		% need to generate matrix of same size as trials group

for i = 1:length(category_group)
	temp = (trials_group == category_group(i));
	temp2 = ones(sum(temp),1)*i;
 	temp_matrix(temp) = temp2;		% assign ranks within temp matrix for each category group
end
trials_group = temp_matrix;		% assigning ranks at end obviates problem of interference between ranks and remaining values

[P_value, ANOVATAB] = anova1(data_cnts, trials_group);

MSgroup = ANOVATAB{2,4};
MSerror = ANOVATAB{3,4};

close
close