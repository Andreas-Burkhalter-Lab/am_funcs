% Test function
function analysis_out = test_seea_func(ephysin, analysis_in)
   analysis_out = analysis_in;
   analysis_out.test = ephysin(1).channels;
end