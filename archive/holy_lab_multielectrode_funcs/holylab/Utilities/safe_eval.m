function eval_out_g=safe_eval(str, vars_to_out_g)
% eval a string w/o worrying about corrupting current namespace
% 
% SYNTAX:
%   eval_out_g=safe_eval(str, vars_to_out_g)
% PRE:
%   str: a string to eval
%   vars_to_out_g: a cell array of variable names to return
% POST:
%   eval_out_g: a struct to hold return values

   if(~iscell(vars_to_out_g))
      vars_to_out_g={vars_to_out_g};
   end
   
   % TODO: check eval_out_g/vars_to_out_g/tttIdx are not in vars_to_out_g

   eval(str);
   
   % hopefully tttIdx is not what user asks to return
   for tttIdx=1:length(vars_to_out_g)
      if(~exist(vars_to_out_g{tttIdx}, 'var'))
         error('no such var exists');
      end
      
      eval_out_g.(vars_to_out_g{tttIdx})=eval(vars_to_out_g{tttIdx});
   end
   