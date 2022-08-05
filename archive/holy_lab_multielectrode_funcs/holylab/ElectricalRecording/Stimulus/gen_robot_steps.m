function gen_robot_steps(nFractions, nMixtures, filename, options)
% usage: 
%   gen_robot_steps(nFractions, nMixtures, filename, options)
% pre: 
%   nFractions: # of fractions a mixture can have
%   nMixtures: # of mixtures, it is also # of available/candidate fractions
%   filename: the output filename. If omitted or empty, output to screen
%   options: a structure which has fields: 
%      LF: the line ending
%      filter_file: the .mat filename that has result of f=mixing_filter(nFractions);
%         if this filename is empty, user will be asked.
%         (Note: (1) the saved variable name must be f;
%                (2) the nFractions used in mixing_filter() must be same as 
%                    in gen_robot_steps() )
%      method: the method name, default is 'Mixing1'. Note: no tab is allowed.
%      sample: the sample name, default is '4mL Sample'. Note: no tab is allowed.
% 
% eg:
%   gen_robot_steps(2, 3)

   if(nargin==4 && isfield(options, 'filter_file'))
      if(isempty(options.filter_file))
         [ttfilename, ttpathname] = uigetfile('*.mat', 'Pick a MAT-file');
         if(isequal(ttfilename,0) | isequal(ttpathname,0))
           disp('User pressed cancel');
           return;
         else
            options.filter_file= fullfile(ttpathname, ttfilename);
         end
      end
      load(options.filter_file); % 
   else
      f=mixing_filter(nFractions);
   end
   
   A=mixing_circular(f, nMixtures);
   mixture_matrix=zeros(nMixtures, nFractions);
   for row=1:nMixtures
      mixture_matrix(row,:)=find(A(row,:)==1);
   end

   if(nargin==2 || isempty(filename))
      fid=1;
   else
      fid=fopen(filename, 'w');
      if(fid<0)
         error(['failed to open file with write permission: ' filename]);
      end
   end

   if(nargin==4 && isfield(options, 'LF'))
      LF=options.LF;
   else
      LF='\n';   
   end

   tHeader='Needle Count\tVariable Count\tMethod Name\tMethod Mode\tTube Number\tDescription\tValid Step\tVar1: Value\tVar1: Name\tVar1: Unit\tVar1: Type\tVar1: Format\t';

   fprintf(fid, tHeader);
   fprintf(fid, LF);
   
   if(nargin==4 && isfield(options, 'method'))
      tMethodName=options.method;
   else
      tMethodName='Mixing1';
   end
   
   if(nargin==4 && isfield(options, 'sample'))
      tSampleName=options.sample;
   else
      tSampleName='4mL Sample';
   end

   for col=1:nFractions
      for row=1:nMixtures
         tTube=row;
         tSamp=mixture_matrix(row, col);
         fprintf(fid, '1\t1\t%s\t0\t%d\t\tT\t%s:%d\tSOURCE\tNO_UNIT\tTArrayTube\tLeft\t', tMethodName, tTube, tSampleName, tSamp);
         fprintf(fid, LF);
      end % for, row
   end % for, col

   if(fid>2)
      fclose(fid);
   end
   
