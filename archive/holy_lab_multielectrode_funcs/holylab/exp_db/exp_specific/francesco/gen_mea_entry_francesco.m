function entry=gen_mea_entry_francesco(analyze_script_filename, options)

  default_options('date', []);
  entry = gen_mea_entry(analyze_script_filename, options);
   
  entry.investigator = {'francesco'};

