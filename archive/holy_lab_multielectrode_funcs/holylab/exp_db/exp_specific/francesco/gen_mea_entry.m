function entry=gen_mea_entry(analyze_script_filename, options)

   default_options('date', []);
   

   [pdir] = fileparts(analyze_script_filename);
   if isempty(pdir)
     pdir = pwd;
   end

   entry.version = 1;
   
   entry.investigator = {getenv('USER')};  % by default, user is investigator

   if(~isempty(options.date))
      entry.start_date = options.date;
   else
      entry.start_date = guess_date(pdir);
      if(isempty(entry.start_date))
         error(['please specify the date for ' analyze_script_filename]);
      end
   end

   % What kind of experiment did you do?
   tags = {'mea','vno'};
   entry.tags = tags;
    
   % Where is your data for this experiment stored?
   entry.data_locations = {pdir};
   
   % now the big job, parse stimulus info
   [tt, filecontent]=load_text_file(analyze_script_filename);
   tExpression= '^(valvelabels\s*=\s*\{.*?\}\;)\s*$'; % lazy mode
   [s, f, t]=regexp(filecontent, tExpression, 'once', 'lineanchors');
   if(isempty(t))
      oldpwd=pwd;
      cd(pdir);
      tt=safe_eval('analyze1', 'valvelabels');
      cd(oldpwd);
   else
      tt=safe_eval(filecontent(t(1):t(2)), 'valvelabels');
   end
   labels=tt.valvelabels;

   stimuli = parse_valvelabels(labels);  % this does the serious work
   entry.stimulus = stimuli;
   
   
function dv=guess_date(pdir)
   dv=[]; % date vector
   
   texpr='Y(\d){2}-?(\d){2}-?(\d){2}';
   [from, to, tokenExtents, match, tokens]=regexp(pdir, texpr, 'once');
   if(~isempty(tokenExtents))
      dv=cellfun(@str2num, tokens);
      dv(1)=2000+dv(1);
      return
   end
   texpr=[filesep '(\d){4}-(\d){2}-(\d){1,2}'];
   [from, to, tokenExtents, match, tokens]=regexp(pdir, texpr, 'once');
   if(~isempty(tokenExtents))
      dv=cellfun(@str2num, tokens);
      return
   end

   
   
   
   
   
   
   
   
   
   
   
   
   

