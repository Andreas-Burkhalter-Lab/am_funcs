function h=show_help(sender, event_args, title, msg)
%  h=show_help(sender, event_args, title, msg)
% USAGE:
%    called directly: show_help(gcf, [], 't', 'm');
%    used in callback: SEE: bind_shortcut()
   h=msgbox(msg, title, 'help', struct('Interpreter', 'tex', 'WindowStyle', 'non-modal'));
   
   
   