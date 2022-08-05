#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include "PointServerImageV.h"
#include "AccumulatorMoments.h"
#include "CriticalNeighborhoodBase.h"
#include "mat.h"
#include "MatlabIO.h"

#include <valgrind/callgrind.h>
#include <time.h>

// Compile with make

using namespace std;
using namespace Eigen;

int main()
{
  // Specify the types we'll be using
  typedef uint8_t ImDataType;
  typedef float   PositionType;
  typedef Eigen::Matrix<PositionType,2,1> CoordP;
  typedef Eigen::Matrix<ImDataType,Eigen::Dynamic,1> CoordV;
  typedef Eigen::Matrix<PositionType,Eigen::Dynamic,1> CoordVP;
  typedef Eigen::Matrix<int,CoordV::RowsAtCompileTime,1> CoordVint;
  typedef PointServerImageV<CoordP,CoordV> PointServerType;
  typedef PointServerType::PointType PointType;
  typedef AccumulatorMoments<PointType> AccumulatorType;

  // The p-value
  PositionType pvalue = 0.001;

  // Load the data
  MatlabIO::Load reader("CN_Image_input.mat");
  mxArray *im = reader.load("im");

  // Prepare all the other information
  const int *pSz = mxGetDimensions(im);
  int n_values = pSz[0];
  Eigen::Vector2i sz,sztmp;
  int i,j;

  sz[0] = pSz[1];
  sz[1] = pSz[2];
  int n_max = sz.maxCoeff();
  int n_points = sz[0]*sz[1];
  
  // Create the PointServer, Accumulator, and PointCollection
  GridNeighbors<CoordP> gn(sz,n_max);
  PointServerType ps(gn,(ImDataType*) mxGetData(im),n_values);
  AccumulatorType acc(n_values);
  CriticalNeighborhoodBase<PointServerType,AccumulatorType> cn(ps,acc,false);
  NeighborhoodHistory history;

  // Initialize the output
  Matrix<PositionType,Dynamic,Dynamic> chisqImage(sz[0],sz[1]);
  Matrix<PositionType,Dynamic,Dynamic> chisqImage1(sz[0],sz[1]);
  Matrix<int,Dynamic,Dynamic> nbrhoodSize(sz[0],sz[1]);
  Matrix<int,Dynamic,Dynamic> nImage(sz[0],sz[1]);
  Matrix<ImDataType,Dynamic,Dynamic> filteredImage(n_values,sz[0]*sz[1]);
  nImage.setZero();
  chisqImage.setZero();
  chisqImage1.setZero();
  nbrhoodSize.setZero();

  // Run 1 point
  sztmp[0] = 79;
  sztmp[1] = 378;
  ps.basePositionG(sztmp);
  CriticalNeighborhoodAlgorithms::flowToPeak(cn,history,pvalue);
  for (i = 0; i < cn.sortOrder().size(); i++) {
    j = cn.sortOrder()[i];
    chisqImage1(j) += cn.chisqNeighborhood()[i];
  }

  // Run points
  vector<PositionType> corr;
  vector<int> orderIndex;
  for (int pixIndex = 0; pixIndex < ps.N(); pixIndex++) {
    if (pixIndex%1000 == 0)
      cout << pixIndex << "..." << flush;
    ps.basePositionIndex(pixIndex);
    int nNbrs = CriticalNeighborhoodAlgorithms::flowToPeak(cn,history,pvalue);
    for (i = 0; i < cn.sortOrder().size(); i++) {
      j = cn.sortOrder()[i];
      nImage(j)++;
      //chisqImage(j) += cn.chisqNeighborhood()[i];
    }
    // Assign only to most correlated points
    CoordVP cortmp;
    CriticalNeighborhoodAlgorithms::computeCorrelations(corr,cortmp,acc,ps,cn.sortOrder());
    orderIndex.resize(cn.sortOrder().size());
    for (i = 0; i < orderIndex.size(); i++)
      orderIndex[i] = i;
    IndirectGreater<vector<PositionType> > comp(corr);
    sort(orderIndex.begin(),orderIndex.end(),comp);
    for (i = 0; i < orderIndex.size()-nNbrs; i++) {
      j = cn.sortOrder()[orderIndex[i]];
      chisqImage(j)++;
    }
    filteredImage.col(pixIndex) = acc.mean().cast<ImDataType>();
    nbrhoodSize(pixIndex) = nNbrs;
  }
  cout << "done." << endl;

  MatlabIO::Save writer("watershed_image_output.mat");
  writer.save(nImage,"n");
  writer.save(chisqImage,"chisq");
  writer.save(chisqImage1,"chisq1");
  writer.save(filteredImage,"filteredImage");
  writer.save(nbrhoodSize,"neighborhoodSize");
}
