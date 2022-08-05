# Note: the following environment variables must be defined:
#   MATLABFUNCDIR
#   MATLABINSTALLDIR
#   EIGENDIR
#   BOOSTDIR
SRCDIR = src
BOOSTSRCDIR = $(BOOSTDIR)/libs/thread/src/pthread
UTILITYDIR = $(MATLABFUNCDIR)/Utilities
IMAGEDIR = $(MATLABFUNCDIR)/images
MATOPSDIR = $(MATLABINSTALLDIR)/bin
CNDEPS = $(SRCDIR)/AccumulatorMoments.h $(SRCDIR)/CriticalNeighborhoodAlgorithms.h $(SRCDIR)/CriticalNeighborhoodBase.h $(SRCDIR)/NeighborhoodStatistics.h $(SRCDIR)/NeighborhoodHistory.h
MATLABDEPS = $(UTILITYDIR)/MatlabIO.h $(UTILITYDIR)/MatlabTraits.h
CFLAGS = -g -O -DMAIN

cn_mclinesearch: cn_mclinesearch.cpp $(SRCDIR)/MonteCarloLineSearch.cpp $(SRCDIR)/T2Direct.h
	mex -O -I$(EIGENDIR) -I$(SRCDIR) -I$(BOOSTDIR) -I$(UTILITYDIR) cn_mclinesearch.cpp $(BOOSTSRCDIR)/thread.cpp $(BOOSTSRCDIR)/once.cpp -lut
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(EIGENDIR) -I$(SRCDIR) -I$(BOOSTDIR) -I$(UTILITYDIR) cn_mclinesearch.cpp $(BOOSTSRCDIR)/thread.cpp $(BOOSTSRCDIR)/once.cpp -lut -o cn_mclinesearch



cn_preordered: cn_preordered.cpp $(SRCDIR)/PointServerPreorderedC.h $(SRCDIR)/AccumulatorMoments.h $(SRCDIR)/CriticalNeighborhoodAlgorithms.h $(MATLABDEPS)
	mex -O -I$(UTILITYDIR) -I$(SRCDIR) cn_preordered.cpp

cn_filter_boxcar: cn_filter_boxcar.cpp $(SRCDIR)/PointServerFilter1d.h $(CNDEPS) $(MATLABDEPS)
	mex -O -I$(UTILITYDIR) -I$(SRCDIR) cn_filter_boxcar.cpp

all: cn_preordered cn_filter_boxcar


# Notes on how to install boost library with the install package from boost.org
# execute the fellowing command in terminal before compiling boost library
#
#     export CFLAGS="$CFLAGS -fPIC"
#     export CXXFLAGS="$CXXFLAGS -fPIC"
#
# These two commands will add -fPIC to gnu compiler.
# After adding this flag & recompiling the boost library, it works. No more complain about -fPIC error.


