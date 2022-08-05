function h = mex_interrupt(str)
  % mex_interrupt: create a message window; dismissal can be tested to stop MEX execution
  h = msgbox(sprintf('Press OK to interrupt execution of MEX file %s. It may take some time to respond.',str),sprintf('Interrupt %s',str),'warn');
end
