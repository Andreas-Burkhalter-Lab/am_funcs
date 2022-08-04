% Scripts for communicating with diesel remotely from Matlab.

%%%%%% line 8 not successfully running 'save' command of testoutput; maybe
%%%%%% create shell script with ssh matlab run command on n5, call script
%%%%%% from HP ssh..... line 12 works

% remotely run a Matlab function
[junk b] = system('plink -pw abcd1234 illya@n5.wustl.edu ssh illya@diesel matlab -nodisplay -r testoutput')
% date =  b(end-62:end-6);

% [junk b] = system('plink -pw abcd1234 illya@n5.wustl.edu ssh -X -t -t illya@diesel rbash /home/illya/Andrew/recordings/run_matlab_outputtest')

% put values from a remote file into the matlab variable 'dosout'
[junk sysout] = system(['plink -pw abcd1234 illya@n5.wustl.edu ssh -X -q -t -t illya@diesel head /home/illya/Andrew/recordings/atest 2>&1 2> deletthis'])

% print the header of the file we created on the remote computer
% dosout = str2num(dosout)

clear junk