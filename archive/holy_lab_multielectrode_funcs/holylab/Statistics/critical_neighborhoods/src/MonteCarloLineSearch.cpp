#include <iostream>
#include <vector>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/thread.hpp>
#include "T2Direct.h"
#include "IndirectCompare.h"
//The following are for debugging, e.g., exporting the data matrix to Matlab
//#include "mat.h"
//#include "MatlabIO.h"

// To compile without Matlab support:
// g++ -g -O -I~/src/eigen/ MonteCarloLineSearch.cpp -lboost_thread-mt -o MonteCarloLineSearch

// To compile with Matlab support
// Because Matlab uses libboost_thread internally, and different
// versions of Matlab use different ABI-incompatible versions of
// libboost_thread, the best option is to statically link:
//    ln -s /usr/lib/libboost_thread-mt.a some-unique-name.a
//    mex -f /usr/local/matlab2011b/bin/matopts.sh -g -O -I~/src/eigen/ -I/home/tim/matlabfunc/Utilities/ -L. MonteCarloLineSearch.cpp -lsome-unique-name -o MonteCarloLineSearch
//
// Another option is to use the same boost version as matlab, discover
// this using
//    locate libboost_thread.so.1
// from the Unix command line. If necessary, download a copy of the
// library and unpack in /usr/local/include. (No other installation is
// needed, even though "threads" usually requires installation,
// because it will use Matlab's installed version of the library.)
// Compilation command on diva (matlab 2011b, using boost 1.44):
// mex -f /usr/local/matlab2011b/bin/matopts.sh -g -O -I/usr/local/include/boost_1_44_0/ -I~/src/eigen/ -I/home/tim/matlabfunc/Utilities/ MonteCarloLineSearch.cpp -o MonteCarloLineSearch

//
// Utility function: take the square of a number
//
template <typename T>
inline T square(const T& value)
{
  return value * value;
}

//
// Given an "inner" set of n points that are the closest n points, and
// another set of "outer" points that are all more distant than the
// inner points, calculate the maximum displacement from the origin (in
// a particular direction) consistent with the inner points remaining
// the closest n points
//
template <typename Derived,typename TempType>
typename Derived::Scalar lineSearchBound(const Derived& xInner,const Derived& ySquaredInner,const Derived& xOuter,const Derived& ySquaredOuter,int s,TempType &tmp)
{
  typename Derived::Scalar alpha;
  int IndexInner,IndexOuter,IndexInnerOld,IndexOuterOld,i;
  const typename Derived::Scalar *p;
  // Choose the inner point based on the flag
  if (s < 0) {
    p = std::max_element(&xInner[0],&xInner[0]+xInner.size());
    IndexInner = p-&xInner[0];
  } else {
    p = std::min_element(&xInner[0],&xInner[0]+xInner.size());
    IndexInner = p-&xInner[0];
  }
  // The outer index starts with the closest point that is in the
  // correct direction
  for (IndexOuter = 0; IndexOuter < xOuter.size(); IndexOuter++)
    if (s*xOuter[IndexOuter] > s*xInner[IndexInner])
      break;
  if (IndexOuter == xOuter.size())
    return s*std::numeric_limits<double>::infinity();
  // Find the cutoff iteratively, always in terms of a chosen pair of
  // points and then checking whether each member of the pair is
  // correctly-chosen
  IndexInnerOld = IndexOuterOld = -1;
  while (IndexInner != IndexInnerOld || IndexOuter != IndexOuterOld) {
    // Store the previous estimates of the constraining point indices
    IndexInnerOld = IndexInner;
    IndexOuterOld = IndexOuter;
    // Compute the midpoint along the line
    alpha = (square(xOuter[IndexOuter]) + ySquaredOuter[IndexOuter] - 
	     square(xInner[IndexInner]) - ySquaredInner[IndexInner])/
      (2*(xOuter[IndexOuter] - xInner[IndexInner]));
    // Compute the squared distance to all inner points, and find the
    // most distant of these
    tmp.leftCols(xInner.size()) = (xInner.array() - alpha).square() + ySquaredInner.array();
    p = std::max_element(&tmp[0],&tmp[0]+xInner.size());
    IndexInner = p-&tmp[0];
    // Compute the squared distance to all outer points, and find the
    // closest of these
    tmp.leftCols(xOuter.size()) = (xOuter.array() - alpha).square() + ySquaredOuter.array();
    p = std::min_element(&tmp[0],&tmp[0]+xOuter.size());
    IndexOuter = p-&tmp[0];
  }
  return alpha;
}

//
// The principal function: calculate T^2 using Monte Carlo methods
//
// The API for this function is deliberately designed to facilitate
// threading: different threads should get different offsets in the T2
// and T2ls arrays, and the RandomEngine should be seeded with
// something like the thread number to ensure that different threads
// generate different random numbers
template <typename T2Type,typename RandomEngineType>
void MonteCarloLineSearch(double *pT2,double *pT2ls,int nsim,const std::vector<int>& nlist,int d,int npoints,T2Type& T2D,RandomEngineType &rng)
{
  // Set up the normal distribution generator
  boost::normal_distribution<> nd(0.0, 1.0);
  boost::variate_generator<RandomEngineType&, 
                           boost::normal_distribution<> > var_nor(rng, nd);

  //
  // Allocate all the storage needed for the algorithm
  //
  int simIndex,nbrhoodSizeIndex,i,jleft,jright;
  double munorm,alphaleft, alpharight,factor;
  // The array holding the simulated data points
  Eigen::MatrixXd X(d,npoints);
  // The data points in sorted order
  Eigen::MatrixXd sX(d,npoints);
  // The square distance
  Eigen::RowVectorXd squaredist(npoints);
  // Sorting the points by distance
  IndirectLess<Eigen::RowVectorXd> comp(squaredist);
  std::vector<int> sortOrder(npoints);
  // The mean
  Eigen::VectorXd mu(d);
  // The covariance
  typename T2Type::CovarianceType C;  // size will be allocated on first use
  // The parallel and perpendicular projections
  Eigen::RowVectorXd proj(npoints);
  Eigen::RowVectorXd perpSquared(npoints);
  Eigen::RowVectorXd xLeft(npoints);
  Eigen::RowVectorXd ySquaredLeft(npoints);
  Eigen::RowVectorXd xRight(npoints);
  Eigen::RowVectorXd ySquaredRight(npoints);
  Eigen::RowVectorXd tmp(npoints);
  // The first coordinate vector, used only in cases where mu is zero
  Eigen::VectorXd e1(d);
  e1.setZero();
  e1(0) = 1;
  
  //
  // Loop over simulations
  //
  for (simIndex = 0; simIndex < nsim; simIndex++) {
    // Generate the data points
    for (i = 0; i < X.size(); i++)
      X(i) = var_nor();
    //MatlabIO::Save writer("mc.mat");
    //writer.save(X,"x");
    // Compute the distance from zero
    squaredist = X.colwise().squaredNorm();
    // Put them in sorted order
    for (i = 0; i < sortOrder.size(); i++)
      sortOrder[i] = i;
    std::sort(sortOrder.begin(),sortOrder.end(),comp);
    for (i = 0; i < sortOrder.size(); i++)
      sX.col(i) = X.col(sortOrder[i]);
    //
    // Loop over neighborhood sizes
    //
    for (nbrhoodSizeIndex = 0; nbrhoodSizeIndex < nlist.size(); nbrhoodSizeIndex++) {
      int n = nlist[nbrhoodSizeIndex];
      if (n < T2D.nMin()) {
	*pT2 = *pT2ls = std::numeric_limits<double>::infinity();
	pT2++;
	pT2ls++;
	continue;
      }
      // Compute T^2
      *pT2 = T2D.calculate(sX.leftCols(n),mu,C);
      // Pick a projection direction (to mimic linesearch)
      munorm = mu.squaredNorm();
      if (munorm == 0)
	mu(0) = 1;
      else {
	munorm = sqrt(munorm);
	mu /= munorm;
      }
      // Compute the projection, and also the perpendicular distance
      // (these are the rho, z components of the coordinates in
      // cylindrical coordinates)
      proj.noalias() = mu.transpose()*sX;
      for (i = 0; i < npoints; i++)
	perpSquared(i) = squaredist(sortOrder[i]) - square(proj(i));
      // Split into center, left, and right groups
      // The center group is just the first n entries
      jleft = jright = 0;
      for (i = n; i < npoints; i++)
	if (proj(i) < 0) {
	  xLeft(jleft) = proj(i);
	  ySquaredLeft(jleft) = perpSquared(i);
	  jleft++;
	} else {
	  xRight(jright) = proj(i);
	  ySquaredRight(jright) = perpSquared(i);
	  jright++;
	}
      // Find the positions along the line, in each direction, for which the
      // distance to the farthest "inside" point is equal to the distance to
      // the closest "outside" point. This defines the span of movements
      // along the line for which the inside points are the closest n
      // points.
      alphaleft = lineSearchBound(proj.leftCols(n),perpSquared.leftCols(n),xLeft.leftCols(jleft),ySquaredLeft.leftCols(jleft),-1,tmp);
      alpharight = lineSearchBound(proj.leftCols(n),perpSquared.leftCols(n),xRight.leftCols(jright),ySquaredRight.leftCols(jright),1,tmp);
      // Compute T^2 for each of these points and pick the largest
      if (munorm > 0) {
	factor = std::max(square(1-alphaleft/munorm),square(1-alpharight/munorm));
	*pT2ls = factor * *pT2;
      } else {
	// The case where mu is zero requires that we evaluate
	// e_1^T C^{-1} e_1
	factor = std::max(square(alphaleft),square(alpharight));
	*pT2ls = factor * T2D.calculate_newmu(e1,C);
      }
      pT2++;
      pT2ls++;
    }
  }
}

//
// Use this if you want to multithread the simulation.
// However, if you're running within Matlab, a better option is to use
// BoostThreadPool.h.
//
template <typename T2Type>
void MonteCarloLineSearchThreaded(double *pT2,double *pT2ls,int nsim,const std::vector<int>& nlist,int d,int npoints,T2Type &T2D,int nthreads)
{
  if (nthreads == 1) {
    boost::mt19937 rng;
    MonteCarloLineSearch(pT2,pT2ls,nsim,nlist,d,npoints,T2D,rng);
  } else {
    int i,offsetl,offsetr,offsetmem;
    // Create separate random number generators for each thread
    boost::mt19937 *rngp = new boost::mt19937[nthreads];
    for (i = 0; i < nthreads; i++)
      rngp[i].seed(i);
    // Create separate T2Direct objects for each thread (because each
    // may require its own internal storage)
    std::vector<T2Type> T2Dv(nthreads,T2D);
    // Create storage for threads
    std::vector<boost::thread *> threadlist;
    boost::thread *threadp;
    // Launch the threads
    for (i = 0; i < nthreads; i++) {
      offsetl = std::floor(i*double(nsim)/nthreads+0.5);
      offsetr = std::floor((i+1)*double(nsim)/nthreads+0.5);
      offsetmem = offsetl*nlist.size();
      threadp = new boost::thread(MonteCarloLineSearch<T2Type,boost::mt19937>,pT2+offsetmem,pT2ls+offsetmem,offsetr-offsetl,nlist,d,npoints,T2Dv[i],rngp[i]);
      threadlist.push_back(threadp);
    }
    // Wait for the threads to finish
    for (i = 0; i < nthreads; i++)
      threadlist[i]->join();
    // Clean up
    for (i = 0; i < nthreads; i++)
      delete threadlist[i];
    delete[] rngp;
  }
}


/*
int main() {
  const int nsim = 1000;
  std::vector<int> nlist;
  T2Direct<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>,Moments::Isotropic> T2D;
  nlist.push_back(2);
  nlist.push_back(3);
  nlist.push_back(5);
  nlist.push_back(10);

  Eigen::MatrixXd T2(nlist.size(),nsim);
  Eigen::MatrixXd T2ls(nlist.size(),nsim);

  // Single-threaded version
  //boost::mt19937 rng;
  //rng.seed(0);
  //MonteCarloLineSearch(&T2(0),&T2ls(0),nsim,nlist,5,100,T2D,rng);

  // Multithreaded version
  MonteCarloLineSearchThreaded(&T2(0),&T2ls(0),nsim,nlist,5,100,T2D,4);

  MatlabIO::Save writer("mc.mat");
  writer.save(T2,"T2");
  writer.save(T2ls,"T2ls");
}
*/
