function autosnip2_AM(opt_autosnip2, opt_snippetfile)
% AUTOSNIP2: a new version of snippetfileauto.m
%%% edited AM 7/13/15 vivid: 
% opt_autosnip2.threshfactor = 1: choose factor by which to multiply threshold 6*mna 


% @pre: 
%    opt_autosnip2:   a structure holding options for this function,
%       whose legal fields are:
%          1. resume=0: if this field is non-zero, resume from a session log file
%                  selected by user using GUI file selector;
%          2. feedback_channels: tell which channels from 0 to 63 are
%                             feedback channels; 
%                             default value depends on hardware (usually 0 for 
%                             chanel, 63 for diesel);
%          3. do_filter=0: if the field is non-zero, do conditional filter before
%                     auto calculate thresholds; 
%          4. compare_hz60=0: if the field is non-zero, compare snippet results w/
%                        and w/o Hz60 filter, and ask user if this filter is
%                        necessory for snipping the group of file;
%          5. snip_range=[-10, 30]: the snip range when cutting snippet for user's
%                                decision on channels  
%          6. time4data_segment=5: decide how long in seconds of data read for calc 
%                               threshhold and psd
%          7. snip_channels='semiauto': legal values are:
%               'semiauto': snip the first file in each group, and popup a window
%                         to have user to select which channels to snip
%               'all': snip all channels
%               a vector of integers: the channels to snip
%          8. flush_valve=0: the flushing valve number 
%          9. quick_snip_time=2: during each stimulus valve is open, how long (in seconds)
%                                snippets are cut before flushing.
%          10. polarity: can be -1, 1, or 0; determines whether program identifies
%                  spikes by looking for a negative deflection, positive
%                  deflection, or either
%          11. extension: what should .merec be replaced with? default is
%                  '.ssnp' (must be provided as string)
%          12. iterative_threshold: if 1, the threshold will be found
%               twice; the first time like normal, then again with all
%               points that exceeded the first threshold pulled out
%               (prevents one set of very large spikes from pulling the
%               threshold so high a second smaller cell is lost; default 0)
%          13.  timestim_calc_vs_intended_tolerance: see timestim
%          14. timestim_wnds: 1st value is wndL and 2nd value is wndR for
%               timestim
% 
%    opt_snippetfile: a structure holding options for snippetfile(),
%                     please see snipoptions() for details. 
%                     @note: its 3 fields (polarity, tofile and outfilename)
%                     are ignored in this (autosnip2()) function.
% @post:
%    3 set of files are created: .env, .ssnp, .vlv; created files are in the same
%    dir as their source .merec files.
% 
% @eg: autosnip2;
%      autosnip2(struct('resume', 1));
% 
% @see: SNIPOPTIONS, SNIPPETFILE

% @history: since 20030827:
%    20071128: removed doall option, because redundant with files option (RCH)
%    20070216: added options to the list of variables saved for resumption (RCH)
%    20061031: added in new options to allow for use snipping ECG info (RCH)
%    20060810: added in  timestim_calc_vs_intended_tolerance option (RCH)
%    20060425: added in .lfpsensor option to allow me to 
%              snippet the right channel automatically 
%              whether or not I have the LFP filter in place (RCH)
%    20060212: 1. added in iterative_threshold option (RCH)
%    20051029: 1. added (in total hack job) opt_autosnip2.doall option
%                 that just processes all .merec files in the directory
%                 organized by name and without asking any questions
%    20040730: 1. allow user to effectively set feedback channels, snip
%                 range, polarity, and extension when calling autosnip2
%              (RCH)
%    20031107: 1. ask user if sort files in each group by time (earliest to latest)
%    20031003: 1. transpose the channel layout, make it consistent w/ merec;
%              2. make tooltip more readable
%              3. don't snip when flush valve is open
%              4. snip at the end of stimulus instead of at the beginning of stimulus
%              5. allow control how long snippets are quickly cut
%    20031002: make the interface(input parameters) consistent
%    20031001: 1. add option so user can snip all channels, or given channels,
%                 or this program pops up a window to have user select channels.
%              2. speed up of snipping the first file in each group
%              3. redesign the channel selection window, and make it clean and convenient:
%                 buttons move to right-bottom; draw less axis labels;
%                 name graph in channel labels; add tooltip to buttons.
% 
%    20030827: 1. add crash resume functionality; in case matlab/autosnip2 crashes,
%                 this will be a big time-saver;
%              2. update comments;

   
   fprintf(['autosnip2: start at ' datestr(now, 0) '\n']);
   
   switch nargin
      case 0
         opt_autosnip2=[]; 
         opt_snippetfile=[];
      case 1
         opt_snippetfile=[];
      case 2
         % do nothing
      otherwise
         error('autosnip2 accept 0~2 arguments');
   end % switch

   % if not present, add this field and set it w/ default value, then
   % later you can use the "field value" as flag;
   % don't use "field is present or not" as flag because it is not good idea.
   % arg 1:
   if(~isfield(opt_autosnip2, 'resume')) % arg 1
      opt_autosnip2.resume=0;
   end
   if(~isfield(opt_autosnip2, 'iterative_threshold'))
     opt_autosnip2.iterative_threshold = 0;
   end
   
   [status, homedir]=system('echo -n $HOME');
   session_dir=[homedir filesep '.autosnip2'];
   if(~isdir(session_dir))
      ttSuc=mkdir('/', session_dir);
      if(~ttSuc)
         error('fail to create dir for session log');
      end
   end

   if(opt_autosnip2.resume)
      % if resume from crash session
      is_resume=1;
      session_file = UIGetFiles('*.mat', ...
                                'please select ONE session log file', ...
                                session_dir ...
                                );
      if(length(session_file)~=1)
         error('you must select ONE session log file to resume');
      end
      session_file=session_file{1};
      load(session_file); % not: load session_file; @see: help load; help save
      step=stepp; % TODO: a hack
      disp(['resume from ' session_file '@roundCur=' num2str(roundCur) ',fileIdx=' num2str(fileIdx) ...
            ',step=' num2str(step) ',nRounds=' num2str(nRounds)]);
      roundLast=roundCur; fileIdxLast=fileIdx; stepLast=step; % backup these variables
      % if, resume from a crash session
   else
      is_resume=0;
      session_file=[session_dir filesep 'session' datestr(now, 30) '.mat'];

      % the default value of arg2 (feedback_channels) depends on hardware, 
      % it will be handled after readheader()
      %
      % NOTE ON ADDING OPTIONS: if any option will be looked for when
      % autosnip2 resumes after a break, you will need to ALSO add its name
      % to the list of things saved (lookfor "var2save" near end of file)
      
      % AM 6/20/15 added check for thresh_factor
      if(~isfield(opt_autosnip2, 'thresh_factor'))
          opt_autosnip2.thresh_factor = 1; % default to using 6*mna
      end
      
      if(isfield(opt_autosnip2, 'feedback_channels'))    % arg 2
          askfeedback_channels=0;
      else
          opt_autosnip2.feedback_channels=0;
          askfeedback_channels=1;
      end
      if(~isfield(opt_autosnip2, 'do_filter'))    % arg 3
         opt_autosnip2.do_filter=0;
      end
      if(~isfield(opt_autosnip2, 'compare_hz60')) % arg 4
         opt_autosnip2.compare_hz60=0;
      end
      if(~isfield(opt_autosnip2, 'snip_range'))   % arg 5
         opt_autosnip2.snip_range=[-10, 30];
         asksnip_range=1;
      else
         asksnip_range=0;
      end
      if(~isfield(opt_autosnip2, 'time4data_segment')) % arg 6
         opt_autosnip2.time4data_segment=5;
      end
      if(~isfield(opt_autosnip2, 'snip_channels'))     % arg 7
         opt_autosnip2.snip_channels='semiauto';
      end
      if(~isfield(opt_autosnip2, 'flush_valve'))       % arg 8
         opt_autosnip2.flush_valve=0;
      end
      if(~isfield(opt_autosnip2, 'quick_snip_time'))   % arg 9
         opt_autosnip2.quick_snip_time=2;
      end
      if(~isfield(opt_autosnip2, 'polarity'))          % arg 10
         opt_autosnip2.polarity=0;
         askpolarity=1;
      else
         askpolarity=0;
      end
      if(isfield(opt_autosnip2, 'extension'))
          askextension=0;
      else
          opt_autosnip2.extension='.ssnp';
          askextension=1;
      end
      if(~isfield(opt_autosnip2, 'timestim_calc_vs_intended_tolerance'))
          timestim_calc_vs_intended_tolerance = 0;
      else
          timestim_calc_vs_intended_tolerance = opt_autosnip2.timestim_calc_vs_intended_tolerance
      end
      if(~isfield(opt_autosnip2, 'use_alt_header_file'))
          use_alt_header_file = 0;
      else
          use_alt_header_file = opt_autosnip2.use_alt_header_file;
      end
      if(~isfield(opt_autosnip2, 'files'))
          ask_files = 1;
      else
          ask_files = 0;
      end
      opt_autosnip2 = default(opt_autosnip2,'no_env',0);
      opt_autosnip2 = default(opt_autosnip2,'no_redo_vlv',0);
      opt_autosnip2 = default(opt_autosnip2,'do_in_time_order',NaN);
      
      nRounds=0; params4snippetfile=[]; % one "round" means one "group" of files
      
      while(1)
         % step 1 of 4: choose file names in UIGetFiles():
           if ask_files
               merecfiles = UIGetFiles('*.merec', ...
                                     'step 1 of 4: please select Merec data files (click OK when done)'...
                                     );
           elseif iscell(opt_autosnip2.files)
               merecfiles = opt_autosnip2.files;
           else
               nFiles = length(opt_autosnip2.files);
               for nthFile = 1:nFiles
                   merecfiles{nthFile} = [pwd '/' char(opt_autosnip2.files(nthFile))];
               end
           end
         if(length(merecfiles)==0)
            % errordlg('please select at lease 1 file', 'autosnip2', 'modal');
            msgAbort='abort autosnip2';
            msgCont='continue selecting files';
            msgDone='done with selecting files';
            usrResponse=questdlg('please select at least 1 file. What next?', ...
                                 'autosnip2', ...
                                 msgAbort, msgCont, msgDone, msgDone);
            switch(usrResponse)
               case msgAbort 
                  error('user aborted autosnip2');
               case msgCont
                  continue;
               case msgDone
                  break;
               otherwise
                  continue; % if user close questdlg instead of click one of the
                            % three buttons
            end % switch
         end % if, user didn't select any file
         
         % ask user if sort files in the group by time
         if isnan(opt_autosnip2.do_in_time_order)
           usrResponse=questdlg('Do you want to sort files in this group by time (earliest to lastest)?',...
             'autosnip2', ...
             'Yes', 'No', 'Yes');
         else
           usrResponse = opt_autosnip2.do_in_time_order;
         end
         if(strcmp(usrResponse,'Yes'))
            merecfiles=sortfilebytime(merecfiles);
         end % if, user asked for sorting file
         
         if use_alt_header_file == 0
             [header,fid] = readheader(merecfiles{1});
         else
             [header,fid_temp] = readheader(use_alt_header_file);
             fid = fopen(merecfiles{1});
         end
         if(fid==-1)
            % errordlg('please select Merec files', 'autosnip2', 'modal');
            msgAbort='abort autosnip2';
            msgCont='continue selecting files';
            msgDone='done with selecting files';
            usrResponse=questdlg('please select Merec files. What next?', ...
                                 'autosnip2', ...
                                 msgAbort, msgCont, msgDone, msgDone);
            switch(usrResponse)
               case msgAbort 
                  error('use abort autosnip2');
               case msgCont
                  continue;
               case msgDone
                  break;
               otherwise
                  continue; % if user close questdlg instead of click one of the
                            % three buttons
            end % switch
         end % if, error when open merec file 
         
         if(should_use_lfs(merecfiles{1}))
            w = readint16lfs(fid, header.numch, [0 header.scanrate*opt_autosnip2.time4data_segment-1], header.headersize);
         else
            w = fread(fid,[header.numch header.scanrate*opt_autosnip2.time4data_segment],'*int16');
         end
         w_orig=w; % backup raw data
         
         bExplicitlyApplyHz60Filter=0;
         if(opt_autosnip2.compare_hz60)
            if(opt_autosnip2.do_filter)
               so = opt_snippetfile;
               so.Fs = header.scanrate;
               so = snipoptions(so);

               % fprintf('temp filtering begin...\n');            
               
               % filter all available channels:
               w = filterint16(so.condfiltb, so.condfilta, w, 1:header.numch);
               
               % fprintf('temp filtering done\n');
            end
            w=double(w);
            h_woHz60=figure('WindowStyle','normal', ...
                            'numbertitle', 'off',...
                            'name', 'psd w/o applying Hz60 filter' ...
                            ); 
            for i=0:63
               subplot(8, 8, i+1);
               ttChIdx=find(header.channels==i);
               if(isempty(ttChIdx))
                  set(gca, 'visible', 'off');
               else
                  psd(w(ttChIdx, :), 32*2048, header.scanrate);
                  xlabel(''); ylabel('');
                  title(num2str(i));
               end
            end
            
            % repeat the same thing as that w/o Hz60:
            w=w_orig;
            so = opt_snippetfile;
            so.Hz60=1;
            if(opt_autosnip2.do_filter)
               so.Fs = header.scanrate;
               so = snipoptions(so);
               
               % filter all available channels:
               w = filterint16(so.condfiltb, so.condfilta, w, 1:header.numch);
            end
            w=double(w);
            h_wHz60=figure('WindowStyle','normal', ...
                           'numbertitle', 'off',...
                           'name', 'psd w/ applying Hz60 filter' ...
                           ); 
            for i=0:63
               subplot(8, 8, i+1);
               ttChIdx=find(header.channels==i);
               if(isempty(ttChIdx))
                  set(gca, 'visible', 'off');
               else
                  psd(w(ttChIdx, :), 32*2048, header.scanrate);
                  xlabel(''); ylabel('');
                  title(num2str(i));
               end
            end
            
            % now ask for user's choice:
            usrResponse=questdlg('apply Hz60 filter?', ...
                                 'autosnip2', ...
                                 'yes', 'no', 'no');
            switch(usrResponse)
               case 'yes'
                  bExplicitlyApplyHz60Filter=1;
               case 'no'
                  bExplicitlyApplyHz60Filter=0;
               otherwise
                  bExplicitlyApplyHz60Filter=0;
                  % if user close questdlg instead of click one of the two buttons.
                  % treat it as 'no'
            end % switch
            
            if(ishandle(h_wHz60))close(h_wHz60); end
            if(ishandle(h_woHz60))close(h_woHz60); end
         end % if, compare_hz60, i.e. compare w/ and w/o Hz60 filter
         
         if(bExplicitlyApplyHz60Filter)
            opt_snippetfile.Hz60=1;
         end
         
         w=w_orig;
         so = opt_snippetfile;
         if(opt_autosnip2.do_filter)
            so.Fs = header.scanrate;
            so = snipoptions(so);
            
            % filter all available channels:
            w = filterint16(so.condfiltb, so.condfilta, w, 1:header.numch);
         end
         w=double(w);
         % @notes: another way: directly copy w from previous calc
         mn = median(w,2);
         wf = w - repmat(mn,1,size(w,2));
         mna = mean(abs(wf),2);
         thresh = [-opt_autosnip2.thresh_factor*6*mna' + mn';...
             opt_autosnip2.thresh_factor*6*mna'+mn']; %%% AM 6/20/15 added scaling by thresh_factor
         if opt_autosnip2.iterative_threshold == 1
           nCh = size(wf,1);
           for iCh = 1:nCh
             ikill = find( abs(wf(iCh,:)) > (6*mna(iCh)) );
             wf_temp = wf(iCh,:);
             wf_temp(ikill) = [];
             mna_new = mean(abs(wf_temp),2);
             thresh(:,iCh) = [-6*mna_new' + mn(iCh); 6*mna_new'+mn(iCh)];
           end
         end
         % thresh = thresh';
         
         if or(asksnip_range==1,askpolarity==1)
             dlgTitle='step 2 of 4: parameters for snipping';
             prompt={'snip range(e.g. -10, 30) :', ...
                     'polarity(one value from -1, 0, 1):' ...
                    };
             defaultValues={[num2str(opt_autosnip2.snip_range(1)) ', ' num2str(opt_autosnip2.snip_range(2))],num2str(opt_autosnip2.polarity)};
             lineNo=1;
         end
         
         ttSnip_channels=opt_autosnip2.snip_channels;
         if(isnumeric(ttSnip_channels))
             if isfield(opt_autosnip2,'lfpsensor')
                 if opt_autosnip2.lfpsensor ==1
                     if length(header.channels)>2
                         lfplikely = 1;
                         if length(header.channels)>3
                             ecglikely = 1;
                         else
                             ecglikely = 0;
                         end
                         save snipout.mat lfplikely
                     end
                     if isequal(header.channels, [55 62 63])
                         ttSnip_channels = 62;
                     end
                     if isequal(header.channels, [55 62 63 54]);
                         ttSnip_channels = 62;
                     end
%                      if isequal(header.channels, [55 63 62]) % totally confused as to how this came to be - in the middle of an experiment.  this is mostly testing code.
%                          ttSnip_channels = 62;
%                      end
                 end
             end                     
            uncheckedChs=setdiff(header.channels, ttSnip_channels);
            if or(asksnip_range==1,askpolarity==1)
                tUserInputs=inputdlg(prompt,dlgTitle,lineNo,defaultValues);
                tSniprange=eval(['[' tUserInputs{1} ']']);
                tPolarity=str2num(tUserInputs{2});
            else
                tSniprange=opt_autosnip2.snip_range;
                tPolarity=opt_autosnip2.polarity;
            end
%            uiwait(msgbox('skipped because user explicitly provided channels to snip','step 3 of 4: choose channels to snip','modal'));
         elseif(ischar(ttSnip_channels) & strcmp(ttSnip_channels, 'all'))
            uncheckedChs=[];
            if or(asksnip_range==1,askpolarity==1)
                tUserInputs=inputdlg(prompt,dlgTitle,lineNo,defaultValues);
                tSniprange=eval(['[' tUserInputs{1} ']']);
                tPolarity=str2num(tUserInputs{2});
            else
                tSniprange=opt_autosnip2.snip_range;
                tPolarity=opt_autosnip2.polarity;
            end
%            uiwait(msgbox('skipped because user explicitly told me to snip all channels','step 3 of 4: choose channels to snip','modal'));
         elseif(ischar(ttSnip_channels) & strcmp(ttSnip_channels, 'semiauto'))
            fprintf(['autosnip2 is snipping' ...
                     ' the channels of the 1st file...\n']);
            
            if(is_robot_stim(header))
               ttTimeRange=[];
            else
               ttTimeRange=parse_stim_seq(header);
            end
            
            if(isempty(ttTimeRange))
               ttScanRange='all';
            else
               % speed up snipping the first file in each group:
               flush_valve=opt_autosnip2.flush_valve;
               quick_snip_time=opt_autosnip2.quick_snip_time;
               if(ttTimeRange(1,1)==flush_valve)
                  ttTimeRange=ttTimeRange(2:end, :);
               else
                  ttTimeRange=ttTimeRange(1:end, :);
               end
               % here I assume the valve before flushing valve is stimulus
               ttTimeRange=ttTimeRange(find(ttTimeRange(:,1)==flush_valve), 2); % skip flush phase
               % during each stimulus event, only cut the last quick_snip_time second snippets
               ttTimeRange=[ttTimeRange-quick_snip_time ttTimeRange]; 
               if(ttTimeRange(1,1)<0)
                  ttTimeRange(1,1)=0; % in case the first flush only open shorter than quick_snip_time
               end
               ttScanRange=ttTimeRange*header.scanrate; % convert time to scan numbers
               % to avoid off-by-1 problem in snippetfile():
               ttScanRange(end, end)=ttScanRange(end, end)-1;
            end
            so.polarity=0;
            so.tofile=0;
            so.use_alt_header_file = use_alt_header_file;
            [tsnip, snip]=snippetfile(merecfiles{1}, ttScanRange, header.channels, thresh, ...
                                      opt_autosnip2.snip_range, so );
            
            % 4 useful variables: tsnip snip so thresh
            
            if(~isfield(opt_autosnip2, 'feedback_channels')) % arg 2
               % if not present, set by default values:
               
               if(strcmp(key2value(header.wholeheader, 'hardware'), 'dumb-box-0@chanel.wustl.edu'))
                  % this is on chanel:
                  opt_autosnip2.feedback_channels=[0 7 56 63]; 
               elseif(strcmp(key2value(header.wholeheader, 'hardware'), 'dumb-box-1@diesel.wustl.edu'))
                  % this is on diesel:
                  opt_autosnip2.feedback_channels=[54 55 62 63];
               else
                  error('unsupport dumb box');
               end
            end
            
            feedback_channels=opt_autosnip2.feedback_channels;
               
            fprintf(['autosnip2 is snipping' ...
                     ' the channels of the 1st file...done\n']);
            
            
            % step 2-3 of 4: choose channels to snip and set parameters for snip;
            % create a figure w/ an onClose event saving unchecked channels to root's
            % userdata property
            set(0, 'tag', '0');
            h=figure('WindowStyle','normal', ...
                     'numbertitle', 'off',...
                     'name', 'step 3 of 4: choose channels to snip (close figure when done)', ...
                     'DeleteFcn', ...
                     'uncheckedBtns=findobj(gcbo, ''style'', ''toggle'', ''value'', [0]); uncheckedChs=get(uncheckedBtns, ''userdata''); tUserdata=get(0, ''userdata''); tUserdata{1}=uncheckedChs; set(0, ''userdata'', tUserdata); ttInt=str2num(get(0, ''tag'')); set(0, ''tag'',num2str(ttInt+1));'...
                     ); % when this figure close, root's tag increase by 1
            
            for i=0:63
               % subplot(8, 8, i+1); % plot row by row
               ttRow=mod(i, 8);
               ttCol=floor(i/8);
               subplot(8, 8, ttRow*8+ttCol+1); % plot column-by-column
               
               ttChIdx=find(header.channels==i);
               if(isempty(ttChIdx))
                  set(gca, 'visible', 'off');
               else
                  if(~isempty(snip{ttChIdx}))
                     plot(opt_autosnip2.snip_range(1):opt_autosnip2.snip_range(2), snip{ttChIdx});
                     axis tight;
                  end
                  if(i>=56)
                     set(gca, 'xtick', opt_autosnip2.snip_range);
                  else
                     set(gca, 'xtick', []);
                  end
                  set(gca, 'ytick', [min(get(gca, 'YLim')) max(get(gca, 'YLim'))]);
                  % if(mod(i, 8)~=0)
                  %    set(gca, 'ytick', []); 
                  % end
                  % title(num2str(i));
                  tLabels=split_label(key2value(header.wholeheader, 'label list'));
                  tTitle=['ch ' num2str(i) ':' tLabels(ttChIdx, :)];
                  title(tTitle);
                  tTitle=['ch#=' num2str(i) ':label=' tLabels(ttChIdx, :)];
                  
                  posAxis=get(gca, 'position');
                  toggleBtnSize=0.3;
                  % left-top:
                  % posToggleBtn=[posAxis(1), posAxis(2)+posAxis(4)*(1-toggleBtnSize),...
                  %                posAxis(3)*toggleBtnSize, posAxis(4)*toggleBtnSize];
                  % right-bottom:
                  posToggleBtn=[posAxis(1)+posAxis(3)*(1-toggleBtnSize), posAxis(2),...
                                posAxis(3)*toggleBtnSize, posAxis(4)*toggleBtnSize];
                  % preset checkbox's status
                  if(isempty(find(feedback_channels==i)))
                     tbChecked=1; % if not in feedback channel list, preset to checked
                  else
                     tbChecked=0; % if in feedback channel list, preset to unchecked
                  end
                  if(isempty(snip{ttChIdx}))
                     tbChecked=0;
                  end
                  
                  uicontrol(gcf, 'units', 'normalized', 'position', posToggleBtn, ...
                            'style', 'toggle', 'value', [tbChecked], 'string', '', ...
                            'tooltipstring', tTitle, ...
                            ... % 'BackgroundColor', [1 .4 .6], ...
                            'userdata', i);
               end
            end
            
            % popup another parameter window:
            hMoreParamInStep2=nonblock_inputdlg(prompt,dlgTitle,lineNo,defaultValues);
            set(hMoreParamInStep2, 'DeleteFcn', ...
                              'ttInt=str2num(get(0, ''tag'')); set(0, ''tag'',num2str(ttInt+1));'...
                              );% when this figure close, root's tag increase by 1
            
            waitfor(0, 'tag', '2'); % wait for both windows close
            
            tRootUserData=get(0, 'userdata'); % get the result from root's userdata property
            %TEH
            if ~iscell(tRootUserData{1})
                tRootUserData{1} = {tRootUserData{1}};
            end
            uncheckedChs=cell2mat(tRootUserData{1});
            uncheckedChs=sort(uncheckedChs);
            
            tUserInputs=tRootUserData{2};
            tSniprange=eval(['[' tUserInputs{1} ']']);
            tPolarity=str2num(tUserInputs{2});
         else
            error('unexpected option value of snip_channels'); 
         end % end of different cases of opt_autosnip2.snip_channels
         
         [channels,chanIndex] = setdiff(header.channels, uncheckedChs);
         thresh = thresh(:,chanIndex);
         
         % step 4 of 4: give more params:
         if or(askextension==1,askfeedback_channels==1)
             dlgTitle='step 4 of 4: please tell me more parameters (click OK when done)';
             prompt={'extension replacement(i.e. what .merec becomes):', ...
                     'feedback channel (usually 0 for chanel, 63 for diesel):', ...
                    };
             defaultValues={opt_autosnip2.extension, num2str(opt_autosnip2.feedback_channels)};          
             lineNo=1;
             answer=inputdlg(prompt,dlgTitle,lineNo,defaultValues);
             tExtension=answer{1};
             tFeedback=str2num(answer{2});
         else
             tExtension=opt_autosnip2.extension;
             tFeedback=opt_autosnip2.feedback_channels;
         end
         
         so = opt_snippetfile;
         so = snipoptions(so);
         so.polarity=tPolarity;
         so.tofile=1;

         nRounds=nRounds+1;
         % unfortunately, following doesn't work:
         % params4snippetfile(nRounds)=struct('merecfiles', merecfiles, ...
         %                                'channels', channels, ...
         %                                'thresh', thresh, ...
         %                                'sniprange', eval(['[' answer{2} ']']), ...
         %                                'extension', answer{1}, ...
         %                                'so', so ...
         %                                );
         params4snippetfile(nRounds).merecfiles=merecfiles;
         params4snippetfile(nRounds).channels=channels;
         params4snippetfile(nRounds).thresh=thresh;
         params4snippetfile(nRounds).sniprange=tSniprange;
         params4snippetfile(nRounds).extension=tExtension;
         params4snippetfile(nRounds).so=so;
         params4snippetfile(nRounds).feedback=tFeedback;
         
         msgCont='continue selecting files';
         msgDone='done with selecting files';
         if ask_files == 1
           usrResponse=questdlg('Done with this round. What next?', ...
             'autosnip2', ...
             msgCont, msgDone, msgCont);
         else
           usrResponse=msgDone;
         end
         switch(usrResponse)
           case msgCont
             continue;
           case msgDone
             break;
           otherwise
             continue; % if user close questdlg instead of click one of the
                       % two buttons
                           
         end % switch

      end % while,

   end % else, not resume from a crash session

   if(is_resume)
      roundStart=roundLast;
      fileIdxStart=fileIdxLast;
   else
      roundStart=1;
      fileIdxStart=1;
   end

   for roundCur=roundStart:nRounds
      merecfiles=params4snippetfile(roundCur).merecfiles;
      channels=params4snippetfile(roundCur).channels;
      thresh=params4snippetfile(roundCur).thresh;
      sniprange=params4snippetfile(roundCur).sniprange;
      extension=params4snippetfile(roundCur).extension;
      so=params4snippetfile(roundCur).so;
      so.interptimes = 1;
      so.interpsnips = 1;
      so.use_alt_header_file = use_alt_header_file;
      feedback=params4snippetfile(roundCur).feedback;
      
      var2save = {'nRounds', 'params4snippetfile', 'roundCur', 'fileIdx', 'stepp','use_alt_header_file'};
        % this is where new variables that will be looked for after a
        % resume, also, need to be added
      for fileIdx = fileIdxStart:length(merecfiles)
         [pathstr,basename,oldExt] = fileparts(merecfiles{fileIdx});
         if strcmp(pathstr,'') % means didn't start with absolute path and should get
             pathstr = pwd;
         end
         pathstr=[pathstr filesep]; % append / at the end of path. // are ok.
         fprintf(['Now processing ' basename ':\n']);
         if ((~is_resume | (is_resume & stepLast==1) ) & (~opt_autosnip2.no_env) )
            fprintf('(1/3) make envelope ... ');
            stepp=1;
            save(session_file, var2save{:});
            [tStatus, tOutput]=system(['calenv -s 100 -o ' pathstr basename ...
                    '.env '  merecfiles{fileIdx}]);
            fprintf('done\n');
            is_resume=0;
         end
         if(~is_resume | (is_resume & stepLast==2) )
            fprintf('(2/3) cut snippets ... \n');
            stepp=2;
            % save session_file nRounds params4snippetfile roundCur fileIdx step;
            save(session_file, var2save{:});
            so.use_alt_header_file = use_alt_header_file;
            snippetfile(merecfiles{fileIdx},'all',channels,thresh, ...
                        sniprange, ...
                        [pathstr basename, extension], ... %
                        so ...
                        );
            fprintf('      cut snippets ... done\n');
            is_resume=0;
         end
         
         if( (~is_resume | (is_resume & stepLast==3) ) & (~opt_autosnip2.no_redo_vlv) )
            fprintf('(3/3) generate valve file ... ');
            stepp=3;
            % save session_file nRounds params4snippetfile roundCur fileIdx step;
            save(session_file, var2save{:});
            tempopt = struct;
            tempopt.feedback = feedback;
            tempopt.timestim_calc_vs_intended_tolerance = timestim_calc_vs_intended_tolerance;
            if isfield(opt_autosnip2,'timestim_wnds')
                tempopt.timestim_wnds = opt_autosnip2.timestim_wnds;
            end
            if isfield(opt_autosnip2,'no_intended_times')
                tempopt.no_intended_times = opt_autosnip2.no_intended_times;
            end
            if isfield(opt_autosnip2,'total_override')
                tempopt.total_override = opt_autosnip2.total_override;
            end
            merec2vlv([pathstr basename '.vlv'], merecfiles{fileIdx}, tempopt ); %
            fprintf('done\n');
            is_resume=0;
         end
         % no need run cleanstim
      end % for, fileIdx
      
      fileIdxStart=1; % in case resume from crash, next round will start w/ 1st file like usual.
      
   end % for, roundCur

   delete(session_file); % when we are here, no crash happened, so delete session log file

   fprintf(['autosnip2: all done at ' datestr(now, 0) '\n']);

   

