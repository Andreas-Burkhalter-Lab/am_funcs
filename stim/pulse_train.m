 %%% send a pulse train at specified rate, voltage, and duration
 
 pulse_rate_hz = 10; 
 total_duration_sec = 60 * 60; 
 volt_on = 5;
 volt_off = 0; 
 volt_on_duration_sec = 0.01; % how long to spend at volt_on; then return to volt_off
 
 %%% create pulse train as uint8
 
 % start daq session
 
 