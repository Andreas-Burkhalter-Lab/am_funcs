function result=is_robot_stim(h)
% is_robot_stim: check if the stimuli were delivered by the robot
% pre:
%    h: a header
% post:
%    result: 1 if the stimuli were delivered by a robot, 0 otherwise.
% 
   stimField=key2value(h.wholeheader, 'stimulus sequence');
   tag='robot tubes:';
   [result, substr]=is_begin_with(stimField, tag);
