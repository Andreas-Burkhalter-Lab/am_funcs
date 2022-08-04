%%% get scope events for otsquarewave, unidirectional acquisition
% updated 4/5/18 on thermaltake
function scopeevents_thisplane = scope_events_otsquarewave(scope_events,iplane)

    go_on = input('Otsquarewave, unidirectional, [1 1 2 2] acquisition - are you sure? (y=yes) ','s');
    if ~strcmp(go_on,'y')
        error('quitting this function')
    end


if iplane == 1
    inds = sort([1:8:height(scope_events), 3:8:height(scope_events)]);
    scopeevents_thisplane = scope_events(inds,:);
elseif iplane == 2
    inds = sort([5:8:height(scope_events), 7:8:height(scope_events)]);
    scopeevents_thisplane = scope_events(inds,:); 
end