#include <iostream>
#include <sstream>
#include <vector>
#include "AccumulatorImagePV.h"
#include "ImageLandmarkGrid.h"
#include "CriticalNeighborhoodBase.h"
#include "NeighborhoodHistory.h"
#include "mat.h"
#include "MatlabIO.h"

#include <time.h>

// Compile with make

// Syntax:
//    cluster_stack_v1 settingsfilename datafilename outputfilename
// See the "Usage line" below for important details

using namespace std;
using namespace Eigen;

// Flow probe point to the peak or to the nearest other probe
// point. This is a lot like
// CriticalNeighborhoodAlgorithms::flowToNeighbor, except we have to
// test to see whether the spatial position is sufficiently changed to
// make it to the neighboring probe point.
template <typename CNType,typename PointType>
CriticalNeighborhoodAlgorithms::FlowResults flowToImLandmark(CNType& cn,NeighborhoodHistory& history,PointType& startingPoint,typename PointType::Scalar minDisplacement,typename CNType::Scalar pvalue,int nMin = 1)
{
  int nNbrs,firstIndex;
  bool isAtMax = false;
  
  history.restart();
  startingPoint = cn.pointServer().basePoint().position();
  while (true) {
    // Compute the critical neighborhood
    nNbrs = cn.inflateToCriterion(pvalue,nMin);
    if (nNbrs < nMin)
      nNbrs = nMin;
    cn.backtrack(nNbrs);
    // Determine whether we have traveled far enough
    if ((cn.pointServer().basePoint().position() - startingPoint).norm() > minDisplacement)
      break;
    // Determine whether we are at the peak
    history.add(cn.sortOrder(),nNbrs);
    isAtMax = history.isAtMax();
    if (isAtMax)
      break;
    // Shift the base point (and perform any other updating)
    cn.pointServer().update(cn.accumulator());
  }

  CriticalNeighborhoodAlgorithms::FlowResults results;
  results.nNbrs = nNbrs;
  results.isAtMax = isAtMax;
  return results;
}

// Utility functions
template <typename T>
void ind2sub(MatrixBase<T>& coord,int index,const MatrixBase<T>& cumsz)
{
  for (int i = coord.size()-1; i >= 0; i--) {
    coord(i) = index/cumsz(i);
    index -= coord(i)*cumsz(i);
  }
}

template <typename T>
int sub2ind(const MatrixBase<T>& coord,const MatrixBase<T>& cumsz)
{
  return (coord.array() * cumsz.array()).sum();
}


int main(int argc,char* argv[])
{
  // Specify the types we'll be using
  //
  // The base type for representing spatial position in "physical"
  // units (e.g., microns). Must be floating-point.
  typedef float PositionType;
  // A vector definining a spatial position.
  // "Matrix" comes from Eigen. See their help page at
  // http://eigen.tuxfamily.org/api/index.html (in particular, notice
  // the link to the Eigen-for-Matlab users page and the short
  // reference).
  // Note: You might get better performance by specifying a fixed
  // spatial dimensionality. I'm leaving this as "Dynamic" for now to
  // facilitate testing on 2-d images.
  typedef Matrix<PositionType,Dynamic,1> CoordP;
  // Your deltaF/F data type. I recommend "single"
  typedef float ImDataType;
  // A vector defining a pixel "value". Dimensionality is "Dynamic"
  // which means you can set it at run-time.
  typedef Matrix<ImDataType,Dynamic,1> CoordV;
  // Now we're getting to the interesting stuff. First, how do we
  // measure distances among pixel values?
  typedef Metric::Euclidean MetricType;      // Euclidean distance
  //typedef Metric::LpDistance<4> MetricType;     // L4-norm based distance
  // How do we find nearby points?
  typedef PointServerImagePV<CoordP,CoordV,MetricType,false> PointServerType;
  // What type of covariance model are we using to model the
  // statistics of pixel values?  Choices are Isotropic, Diagonal, or
  // Full. A good default is Isotropic, but the others are worth
  // considering.
  typedef AccumulatorImagePV<CoordP,CoordV,Moments::Isotropic> AccumulatorType;
  // For now just use the base CriticalNeighborhood class, but we
  // could customize
  typedef CriticalNeighborhoodBase<PointServerType,AccumulatorType> MyCNType;


  // Check command-line arguments
  if (argc < 4) {
    cout << "Usage:\n" << argv[0]
	 << " settingsfilename datafilename outputfilename\n"
	 << "where all files are .mat files." << endl;
    cout << "Input files must have the following variables (may be outdated, check code for accuracy):\n"
	 << "  settings file contains: see code (pvalue, pixel_skip, etc.)\n"
	 << "  data file contains: data" << endl;
    cout << "Important note: the first dimension of data must be the 'value' dimension.\nFor example, an RGB image would be 3-by-height-by-width in size." << endl;
    cout << "Alternative form:\n" << argv[0]
	 << " settingsfilename datafilename outputfilename peaknbrlistfilename\n"
	 << "where the fourth file will hold information about the neighborhood for each landmark flowed to peak" << endl;
    exit(1);
  }

  // Load parameters
  //
  MatlabIO::Load parameterReader(argv[1]);
  // The p-value
  PositionType pvalue;
  parameterReader.loadScalar("pvalue",pvalue);

  // The displacement of adjacent pixels, in physical
  // units. Recommended: use microns, e.g., [0.71 0.71 5] or whatever
  // your settings are.
  Map<Matrix<PositionType,Dynamic,Dynamic> > pixelSkip(NULL,0,0);
  parameterReader.load("pixel_skip",pixelSkip);
  // Make sure pixelSkip is a row vector
  if (pixelSkip.rows() > 1) {
    cout << "pixelSkip must be a row vector." << endl;
    exit(1);
  }

  // The radius, in physical units, of the maximum neighborhood
  // size. You can safely NOT specify this parameter, and just let it
  // choose the default, unless you find that you are running out of
  // memory.  (At a minimum, this had better be at least the size of a
  // glomerulus.)
  PositionType rmax;
  if (parameterReader.hasVariable("rmax"))
    parameterReader.loadScalar("rmax",rmax);
  else
    rmax = -1;  // will calculate based on image size (to include whole image)

  // Set the spatial distance between landmarks (in physical units).
  // This parameter effectively determines how many "probe point"
  // pixels you are going to analyze.  It's analogous to choosing the
  // number of landmarks and doing landmark flow.
  PositionType probePointSkip;
  parameterReader.loadScalar("probe_skip",probePointSkip);

  // The "value coefficient," a factor converting differences in value
  // to a difference in position. The results are likely to be
  // sensitive to this parameter. There are methods for "learning"
  // this automatically (already available in the code tree), and we
  // can implement them fairly easily if you like. But a good start is
  // to set this by hand (and it might even be the best choice).
  //   Be aware that running speed will also be heavily dependent upon
  // this setting, because neighbors are checked in order of physical
  // proximity.  A very high value coefficient therefore means that
  // spatial distance is somewhat irrelevant, and so many pixels will
  // need to be examined before one can be confident that the closest
  // neighbors have been found.  Conversely, a small value coefficient
  // means that spatial displacement is predictive of total distance,
  // and only a few "extra" pixels will need to be checked.
  //   Units: if physical displacement is measured in microns, and
  // intensity is measured fractionally (e.g., dfof = 0.05 for a 5%
  // change in intensity), then the units are microns/dfof.  If your
  // glomeruli are ~10 microns, and 5% is a typical dfof, then setting
  // valueCoefficient = 10/0.05 might be a good start.
  PositionType valueCoefficient;
  parameterReader.loadScalar("value_coefficient",valueCoefficient);

  // Optionally, choose to flow only probe points that have at least
  // one of their coordinates bigger than some threshold. Think "don't
  // worry about unresponsive pixels" here.
  PositionType valueThreshold;
  if (parameterReader.hasVariable("value_threshold"))
    parameterReader.loadScalar("value_threshold",valueThreshold);
  else
    valueThreshold = 0;

  // Set the minimimum number of points so that any "clump" of this
  // size (or larger) is likely to be stable under bootstrap
  // resampling, i.e., not have all points in the clump be omitted
  // from one of the resampled data sets (assuming 1/pvalue
  // resamplings).
  // This is probably mostly a "nicety" and I do not expect it to be a
  // crucial parameter.
  int nMin = ceil(-log(pvalue));

  // Load data
  //
  MatlabIO::Load dataReader(argv[2]);
  mxArray *data = dataReader.load("data");

  // Extract dimensional information
  //
  const int *pSz = mxGetDimensions(data);
  int n_values = pSz[0];
  int n_dimensions = mxGetNumberOfDimensions(data)-1; // -1 cuz value dimension
  int i,j;
  Matrix<int,1,Dynamic> sz(n_dimensions);
  for (i = 0; i < n_dimensions; i++)
    sz[i] = pSz[i+1];
  int n_points = sz.prod();
  // Default value of rmax is the diagonal size of the entire image,
  // meaning that all points can be candidate neighbors (in this case,
  // this parameter has no impact on the answer). The main reason to
  // choose a smaller value is if you have memory limitations.
  if (rmax < 0)
    rmax = sqrt((sz.cast<PositionType>().array() * pixelSkip.array()).square().sum());


  // Display parameter settings for the user
  //
  cout << "pvalue: " << pvalue << endl;
  cout << "pixel_skip: " << pixelSkip << endl;
  cout << "rmax: " << rmax << endl;
  cout << "probe_skip: " << probePointSkip << endl;
  cout << "value_coefficient: " << valueCoefficient << endl;
  if (valueThreshold > 0)
    cout << "value_threshold : " << valueThreshold << endl;

  // Create the image data PointServer
  //
  // First set up a matrix that encodes the geometry: a conversion
  // factor from "pixel" (or "grid") coordinates (integers) to
  // "physical" coordinates (e.g., microns)
  Matrix<PositionType,Dynamic,Dynamic> PixelCoordsToPhysicalCoords(n_dimensions,n_dimensions);
  PixelCoordsToPhysicalCoords.diagonal() = pixelSkip.transpose();
  // Use this matrix to specify which pixels are neighbors of
  // increasing distance from the center
  GridNeighbors<CoordP> gn(sz,rmax,PixelCoordsToPhysicalCoords);
  // Finally, set up the point server for image points
  MetricType metric;
  PointServerType ps(gn,(ImDataType*) mxGetData(data),n_values,metric);
  ps.valueCoefficient(valueCoefficient);

  // Prepare the timer
  clock_t t0,t1,ttmp;

  // Create the landmark "PointServer"
  //
  cout << "Assigning pixels to landmarks..." << flush;
  t0 = clock();  // time how long this takes
  ImageLandmarkGrid<CoordP,CoordV> lmps(gn,probePointSkip,(ImDataType*) mxGetData(data),n_values,valueCoefficient,metric);
  t1 = clock();
  int dt0 = (int)((t1-t0)/CLOCKS_PER_SEC);
  cout << "done." << endl;
  cout << "Landmark pixel separation: " << lmps.landmarkSeparationG().transpose() << endl;

  // Create the output file and save info about landmarks
  //
  MatlabIO::Save writer(argv[3]);
  writer.save(sz,"stacksz");
  writer.save(lmps.pixelLandmarkAssignment(),"landmark_assignment");
  writer.save(lmps.landmarkPixelIndex(),"landmark_pixel_index");
  writer.save(lmps.landmarkSeparationG(),"probe_point_skip");

  // Create the Accumulator
  //
  AccumulatorType acc(n_dimensions,n_values);

  // Create the CriticalNeighborhood object (a "package" of everything else)
  //
  MyCNType cn(ps,acc,false);

  // Calculate the minimum distance that each probe point needs to
  // flow before becoming closer to another probe point's starting
  // point. Allow a little "extra" for safety.
  PositionType probePointMinDisplacement = 1.1 * probePointSkip/2 * sqrt(PositionType(n_dimensions));

  // Allocate the storage for point-flow history
  NeighborhoodHistory history;
  CriticalNeighborhoodAlgorithms::FlowResults flowResults;
  CoordP startingPoint(n_dimensions);

  // Allocate the output
  vector<int> closestNeighborIndex;
  vector<int> nNbrs;
  vector<bool> isAtPeak;
  vector<bool> isFlowed;

  // Open the peak info file, if desired
  bool savingPeakNbrInfo = false;
  char* peakname = NULL;
  if (argc > 4) {
    peakname = argv[4];
    savingPeakNbrInfo = true;
  }
  MatlabIO::Save peakwriter(peakname);
  stringstream ss;


  // Flow each probe point to the local peak, or until it has moved
  // more than half way to another point on the landmark grid
  //
  t0 = clock();
  ttmp = t0;
  cout << "Initial flow of landmarks (" << lmps.nLandmarks() << " total): " << flush;
  for (i = 0; i < lmps.nLandmarks(); i++) {
    // Every 5s, notify the user about progress
    t1 = clock();
    if (t1-ttmp > 5*CLOCKS_PER_SEC) {
      cout << i << "..." << flush;
      ttmp = t1;
    }
    // Set the base point to be one of the landmarks
    ps.basePointIndex(lmps.landmarkPixelIndex(i));
    if ((ps.basePoint().value().array().abs() >= valueThreshold).any()) {
      isFlowed.push_back(true);
      // Flow it to the peak or to the nearest neighboring landmark
      flowResults = flowToImLandmark(cn,history,startingPoint,probePointMinDisplacement,pvalue,nMin);
      // Store the results
      closestNeighborIndex.push_back(cn.sortOrder().front());
      nNbrs.push_back(flowResults.nNbrs);
      isAtPeak.push_back(flowResults.isAtMax);
      if (flowResults.isAtMax && savingPeakNbrInfo) {
	ss.clear();
	ss.str("");
	ss << "lm" << i;
	peakwriter.save(cn.sortOrder(),ss.str().c_str());
      }
    } else {
      isFlowed.push_back(false);
      closestNeighborIndex.push_back(lmps.landmarkPixelIndex(i));
      nNbrs.push_back(-1);
      isAtPeak.push_back(false);
    }
  }
  t1 = clock();
  int dt1 = (int)((t1-t0)/CLOCKS_PER_SEC);
  cout << "done." << endl;
  if (valueThreshold != 0) {
    // Report the total number of points flowed
    i = 0;
    for (vector<bool>::iterator it = isFlowed.begin(); it < isFlowed.end(); it++)
      i += *it;
    cout << i << " points flowed (the rest did not cross threshold)." << endl;
  }

  // Save the output at this stage
  if (valueThreshold != 0)
    writer.save(isFlowed,"is_flowed");
  writer.save(closestNeighborIndex,"closest_neighbor_index0");
  writer.save(nNbrs,"n_neighbors0");
  writer.save(isAtPeak,"is_at_peak0");


  // Create the "map" of probe points.
  //
  // Later we'll use "map flow" to determine stable points
  std::vector<int> map0(closestNeighborIndex.size());
  for (i = 0; i < map0.size(); i++)
    map0[i] = lmps.pixelLandmarkAssignment(closestNeighborIndex[i]);

  writer.save(map0,"map0_0");



  // For each point that appears to be flowing to itself, "reflow" to
  // the true peak (to make certain it really is stable)
  //
  cout << "Starting reflow phase:" << endl;
  std::vector<int> map;
  t0 = clock();
  int nReflowed = 1;
  while (nReflowed) {
    nReflowed = 0;
    // Flow the map
    map = map0;
    CriticalNeighborhoodAlgorithms::flowMap(map,nNbrs);
    // Check each point
    for (i = 0; i < map.size(); i++) {
      if (isFlowed[i] && map[i] == i && !isAtPeak[i]) {
	ps.basePointIndex(lmps.landmarkPixelIndex(i));
	nNbrs[i] = CriticalNeighborhoodAlgorithms::flowToPeak(cn,history,pvalue,nMin);      
	closestNeighborIndex[i] = cn.sortOrder().front();
	map0[i] = lmps.pixelLandmarkAssignment(closestNeighborIndex[i]);
	isAtPeak[i] = true;
	if (savingPeakNbrInfo) {
	  ss.clear();
	  ss.str("");
	  ss << "lm" << i;
	  peakwriter.save(cn.sortOrder(),ss.str().c_str());
	}
	nReflowed++;
      }
    }
    cout << "Reflowed " << nReflowed << " points." << endl;
  }
  t1 = clock();
  int dt2 = (int)((t1-t0)/CLOCKS_PER_SEC);

  writer.save(closestNeighborIndex,"closest_neighbor_index");
  writer.save(nNbrs,"n_neighbors");
  writer.save(isAtPeak,"is_at_peak");
  writer.save(map0,"map0");
  writer.save(map,"map");

  cout << "Landmarking all pixels, total time: " << dt0 << "s" << endl;
  cout << "Initial flow, total time: " << dt1 << "s" << endl;
  cout << "Reflow, total time: " << dt2 << "s" << endl;
}

