warning(['If this fails, edit your matops.sh file to change the fortran' ...
	 ' build from g95 to gfortran']);
mex reorth_mex.c reorth.f -lblas -output reorth
mex tqlb_mex.c tqlb.f -lblas -llapack -output tqlb
mex bdsqr_mex.c dbdqr.f -lblas -llapack -output bdsqr
