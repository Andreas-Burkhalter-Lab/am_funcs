UTILITYDIR = $(MATLABFUNCDIR)/Utilities
MATOPSDIR = $(MATLABINSTALLDIR)/bin
CFLAGS = -g -O -DMAIN -I$(UTILITYDIR) -I/usr/include/valgrind

iminterp: iminterp.cpp $(UTILITYDIR)/image_utilities.cpp image_utilities.h
	mex iminterp.cpp $(UTILITYDIR)/image_utilities.cpp

#imrestrict_mex: imrestrict_mex.cpp image_utilities.cpp image_utilities.h imiterators.cxx
#	mex imrestrict_mex.cpp image_utilities.cpp

imqinterp_st: imqinterp_st.cpp image_utilities.cpp image_utilities.h $(UTILITYDIR)/mat_variables.c $(UTILITYDIR)/mat_variables.h
	mex -I$(UTILITYDIR) imqinterp_st.cpp image_utilities.cpp $(UTILITYDIR)/timer_g.c
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) imqinterp_st.cpp image_utilities.cpp $(UTILITYDIR)/mat_variables.c $(UTILITYDIR)/timer_g.c -o imqinterp_st


imqinterp: imqinterp_mt.cpp image_utilities.cpp image_utilities.h $(UTILITYDIR)/mat_variables.c $(UTILITYDIR)/mat_variables.h
	mex -I$(UTILITYDIR) imqinterp_mt.cpp image_utilities.cpp $(UTILITYDIR)/timer_g.c -output imqinterp
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) imqinterp_mt.cpp image_utilities.cpp $(UTILITYDIR)/mat_variables.c $(UTILITYDIR)/timer_g.c -o imqinterp
