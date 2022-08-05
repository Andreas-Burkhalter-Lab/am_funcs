% add my exp to the db. Run this on envy.wustl.edu

cd /mnt/jason_behavior_e/jason/exp_data_4219to4624.full

rand('twister',5489); % make debug easier

to_validate=false;

entry=gen_xdb_entry_behavior({'20040407', '20040408', '20040414', '20040415'}, struct('phase', 'a'));
if(to_validate)
   if(~validate_xdb_entry_var(entry))
      error('pls check your xdb entry generation function');
   end
end

save_xdb_entry(entry);


entry=gen_xdb_entry_behavior({'20040421', '20040422'}, struct('phase', 'b'));
save_xdb_entry(entry);

entry=gen_xdb_entry_behavior({'20040505', '20040506', '20040512', '20040513'}, struct('phase', 'd'));
save_xdb_entry(entry);

entry=gen_xdb_entry_behavior({'20040519', '20040520'}, struct('phase', 'e'));
save_xdb_entry(entry);

entry=gen_xdb_entry_behavior({'20040623', '20040624'}, struct('phase', 'f'));
save_xdb_entry(entry);


% now test search functions:

% a more flexible way:
result=search_xdb(@query_callback_behavior) % some args can be passed too.

result{1}.data_locations{:}
result{1}.comment

% a more convenient way:
conditions=struct(...
   'investigator', {{@include, 'jason'}}, ...
   'start_date', {{@after, datenum('20040501', 'yyyymmdd'), @before, datenum('20040601', 'yyyymmdd')}} ... 
   );

result2=search_xdb({@query_callback_and, conditions})
for r=result2
   r{1}.comment
end
 