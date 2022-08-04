%%%%% initialize values for the stim parameter being tested
% updated 18/12/13 on thermaltake 

%%%%%%%%%% for experiments before 1 march 2018 with old stimpars format
% exp_date = datetime(pars.starttime);
% if exp_date < datetime('1-March-2018 00:00:00') %%% if experiment split sftf and orient into separate runs, using old format
%     if pars.nAngles < 2 % if sf and tf rather than orient were tested
%         pars.angle_fixed = pars.angle_range(1); 
%         pars.sf_vals = unique(par_sets.sfreq(par_sets.sfreq~=pars.sf_fixed))';
%         pars.tf_vals = unique(par_sets.tfreq(par_sets.tfreq~=pars.tf_fixed))';
%     else
%         pars.angle_vals = unique(par_sets.Angle)'; % no angle_fixed value required
%         if ~isfield(pars,'tf_fixed')
%             pars.tf_fixed = pars.tf_minmax(1);
%         end
%     end
% end
     
switch thispar % check if the parameter was tested for
    case 'sf'
        if isfield(stimpars,'n_sfs') && stimpars.n_sfs > 1
            do_analysis = 1;
            stimparvals = stimpars.sf_vals';
            nparvals = stimpars.n_sfs;
            par_sets_varname = 'sfreq';
            otherparsfixed = stimpar_sets.tfreq == stimpars.tf_fixed & stimpar_sets.orient == stimpars.orient_fixed & rounded_diams == stimpars.diam_minmax(1);
        end
    case 'tf'
        if isfield(stimpars,'n_tfs') && stimpars.n_tfs > 1
            do_analysis = 1;
            stimparvals = stimpars.tf_vals';
            nparvals = stimpars.n_tfs;
            par_sets_varname = 'tfreq';
            otherparsfixed = stimpar_sets.sfreq == stimpars.sf_fixed & stimpar_sets.orient == stimpars.orient_fixed & rounded_diams == stimpars.diam_minmax(1);
        end
    case 'orient'
        if isfield(stimpars,'n_orients') && stimpars.n_orients > 1
            do_analysis = 1;
            stimparvals = stimpars.orient_vals';
            nparvals = stimpars.n_orients;
            par_sets_varname = 'orient';
            otherparsfixed = stimpar_sets.sfreq == stimpars.sf_fixed & stimpar_sets.tfreq == stimpars.tf_fixed & rounded_diams == stimpars.diam_minmax(1);
        end
    case 'diam'
        if isfield(stimpars,'n_diams') && stimpars.n_diams > 1
            do_analysis = 1;
            stimparvals = stimpars.diam_vals';
            nparvals = stimpars.n_diams;
            par_sets_varname = 'diam';
            otherparsfixed = stimpar_sets.sfreq == stimpars.sf_fixed & stimpar_sets.tfreq == stimpars.tf_fixed & stimpar_sets.orient == stimpars.orient_fixed;
        end
end