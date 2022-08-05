UTILITYDIR = $(MATLABFUNCDIR)/Utilities
MATOPSDIR = $(MATLABINSTALLDIR)/bin
CFLAGS = -g
CNDEPENDENCIES = ../AccumulatorMoments.h ../CriticalNeighborhoodAlgorithms.h ../CriticalNeighborhoodBase.h ../NeighborhoodStatistics.h ../NeighborhoodHistory.h
OUTLIERDEPENDENCIES = ../OutlierStatistics.h ../Metric.h
MATLABIODEPENDENCIES = $(UTILITYDIR)/MatlabIO.h $(UTILITYDIR)/MatlabTraits.h

PSFilter: PSFilter.cpp ../PointServerFilter1d.h $(OUTLIERDEPENDENCIES) $(UTILITYDIR)/MatlabIO.h $(UTILITYDIR)/MatlabTraits.h
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) -I.. PSFilter.cpp -o PSFilter

CN_PSFilter: CN_PSFilter.cpp ../PointServerFilter1d.h $(CNDEPENDENCIES) $(OUTLIERDEPENDENCIES) $(MATLABIODEPENDENCIES)
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) -I.. CN_PSFilter.cpp -o CN_PSFilter

PSDistanceDirect: PSDistanceDirect.cpp ../PointServerDistanceDirectC.h ../Metric.h ../OutlierStatistics.h $(UTILITYDIR)/MatlabIO.h $(UTILITYDIR)/MatlabTraits.h
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) -I.. PSDistanceDirect.cpp -o PSDistanceDirect

CN_PSDistanceDirect: CN_PSDistanceDirect.cpp ../PointServerDistanceDirectC.h ../AccumulatorMoments.h ../CriticalNeighborhoodAlgorithms.h ../CriticalNeighborhoodBase.h ../Metric.h $(UTILITYDIR)/MatlabIO.h $(UTILITYDIR)/MatlabTraits.h
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) -I.. CN_PSDistanceDirect.cpp -o CN_PSDistanceDirect

CN_PSPreordered: CN_PSPreordered.cpp ../PointServerPreorderedC.h ../AccumulatorMoments.h ../CriticalNeighborhoodAlgorithms.h ../CriticalNeighborhoodBase.h ../Metric.h $(UTILITYDIR)/MatlabIO.h $(UTILITYDIR)/MatlabTraits.h
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) -I.. CN_PSPreordered.cpp -o CN_PSPreordered

PSGrid: PSGrid.cpp ../PointServerGrid.h ../Metric.h $(UTILITYDIR)/MatlabIO.h $(UTILITYDIR)/MatlabTraits.h
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) -I.. PSGrid.cpp -o PSGrid

PSImage: PSImage.cpp ../PointServerImagePV.h ../PointServerGrid.h ../Metric.h $(UTILITYDIR)/MatlabIO.h $(UTILITYDIR)/MatlabTraits.h
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) -I.. PSImage.cpp -o PSImage

PSImageV: PSImageV.cpp ../PointServerImageV.h ../PointServerGrid.h ../Metric.h $(UTILITYDIR)/MatlabIO.h $(UTILITYDIR)/MatlabTraits.h
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) -I.. PSImageV.cpp -o PSImageV

CN_PSImage: CN_PSImage.cpp ../PointServerImagePV.h ../PointServerGrid.h ../AccumulatorImagePV.h ../AccumulatorMoments.h ../CriticalNeighborhoodAlgorithms.h ../CriticalNeighborhoodBase.h ../Metric.h $(UTILITYDIR)/MatlabIO.h $(UTILITYDIR)/MatlabTraits.h
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) -I.. CN_PSImage.cpp -o CN_PSImage

flow_image: flow_image.cpp ../PointServerImagePV.h ../PointServerGrid.h ../AccumulatorImagePV.h ../AccumulatorMoments.h ../CriticalNeighborhoodAlgorithms.h ../CriticalNeighborhoodBase.h ../Metric.h ../NeighborhoodHistory.h $(UTILITYDIR)/MatlabIO.h $(UTILITYDIR)/MatlabTraits.h
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) -I.. flow_image.cpp -o flow_image

clusterAllFlow: clusterAllFlow.cpp ../PointServerDistanceDirectC.h ../AccumulatorMoments.h ../CriticalNeighborhoodBase.h ../CriticalNeighborhoodAlgorithms.h ../Metric.h ../NeighborhoodHistory.h $(UTILITYDIR)/MatlabIO.h $(UTILITYDIR)/MatlabTraits.h
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) -I.. clusterAllFlow.cpp -o clusterAllFlow

cluster_image: cluster_image.cpp ../PointServerDistanceDirectC.h ../AccumulatorMoments.h ../CriticalNeighborhoodAlgorithms.h ../CriticalNeighborhoodBase.h ../Metric.h ../NeighborhoodHistory.h $(UTILITYDIR)/MatlabIO.h $(UTILITYDIR)/MatlabTraits.h
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) -I.. cluster_image.cpp -o cluster_image

watershed_imageV: watershed_imageV.cpp ../PointServerImageV.h ../AccumulatorMoments.h ../CriticalNeighborhoodAlgorithms.h ../CriticalNeighborhoodBase.h ../Metric.h ../NeighborhoodHistory.h $(UTILITYDIR)/MatlabIO.h $(UTILITYDIR)/MatlabTraits.h
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) -I.. watershed_imageV.cpp -o watershed_imageV

cluster_stack_v1: cluster_stack_v1.cpp ../PointServerImagePV.h ../PointServerGrid.h ../AccumulatorImagePV.h ../AccumulatorMoments.h ../CriticalNeighborhoodAlgorithms.h ../CriticalNeighborhoodBase.h ../Metric.h ../NeighborhoodHistory.h $(UTILITYDIR)/MatlabIO.h $(UTILITYDIR)/MatlabTraits.h
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) -I.. cluster_stack_v1.cpp -o cluster_stack_v1

all: CN_PSDistanceDirect CN_PSPreordered PSGrid PSImage PSImageV CN_PSImage flow_image clusterAllFlow cluster_image watershed_imageV cluster_stack_v1
