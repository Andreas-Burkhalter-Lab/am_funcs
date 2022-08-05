#include "landmarked_neighbors.cpp"
#include "msams_converge_core.h"
#include <pthread.h>
#include <string.h>
#include <assert.h>

// Copyright 2006-2008 by Timothy E. Holy

// This is the function that is called by pthread
// It just casts the void* to an appropriate pointer type
template <class Tdata,class Tint>
void *msams_core_thread(void *p)
{
  msams_core_work((thread_data_type<Tdata,Tint>*) p);
}

// This does housekeeping and launches the threads.  Question: what if
// user wants to start from y that are not on the landmarks?  May need
// some new fields in the option structure, or change the decision re
// lmIndex.
template <class Tdata,class Tint>
int msams_core(const Tdata *x,int d,int N,landmarkStruct<Tdata,Tint> &lminfo,int q,outputStruct<Tdata,Tint> &out,const optionStruct &ops)
{
  //Tdata *ytmp;
  int i,j,thisIndex;
  thread_data_type<Tdata,Tint> td;
  pthread_t *thread;

  // Pack arguments into thread data type
  td.x = x;
  td.d = d;
  td.N = N;
  td.lminfo = &lminfo;
  td.q = q;
  td.out = &out;
  td.ops = &ops;

  // Launch threads
  if (q == 1 || ops.n_threads == 1) {
    // single thread
    td.probePtIndex_start = 0;
    td.probePtIndex_end = q;
    msams_core_work(&td);
  }
  else {
    // multiple threads
    thread_data_type<Tdata,Tint> *tdcopies = new thread_data_type<Tdata,Tint>[ops.n_threads];
    thread = new pthread_t[ops.n_threads];
    for (i = 0; i < ops.n_threads; i++) {
      // Copy the info
      memcpy(&(tdcopies[i]),&td,sizeof(thread_data_type<Tdata,Tint>));
      // Customize the ranges
      tdcopies[i].probePtIndex_start = int(q * (float(i)/ops.n_threads));
      tdcopies[i].probePtIndex_end = int(q * (float(i+1)/ops.n_threads));
      if (pthread_create(thread+i,NULL,msams_core_thread<Tdata,Tint>,(void*) &(tdcopies[i])) != 0)
	return 0;  // error creating thread
    }

    // Wait for threads to return
    for (i = 0; i < ops.n_threads; i++)
      pthread_join(thread[i],NULL);
    delete[] tdcopies;
    delete[] thread;
  }
  return 1;   // success
}

// This does the real work
template <class Tdata,class Tint>
void msams_core_work(thread_data_type<Tdata,Tint> *td)
{
  int probeIndex,d,n_landmarks,i,iter;
  int closestIndex,closestIndexOld,closestDataIndex;
  Tdata *thisY, *thisYend;
  bool isdone,iscycling;
  Tdata *ytraj, *ytrajtmp, *R2traj;
  Tint *ntraj, *xnbrI;
  bool traj_is_tmp;
  int n;
  Tdata R2;

  Tdata tol = td->ops->convergence_thresh;
  

  /*
  // Thread debugging
  for (i = td->probePtIndex_start; i < td->probePtIndex_end; i++)
    td->out->closestDataIndex[i] = (Tint) pthread_self();
  */
  
  // Unpack commonly-used information from the structure created to
  // support multithreading
  d = td->d;
  n_landmarks = td->lminfo->n_landmarks;

  // Allocate temporary storage
  vector<Tdata> x_cum(d);
  vector<Tdata> deltax_cum(d);
  vector<Tdata> x_backtrack(d);
  vector<Tdata> sumsq_cum(d);
  landmarked_neighbors<Tdata,Tint> lm_nbrs;
  lm_nbrs.allocate(n_landmarks,td->N);


  if (td->out->ytraj == NULL) {
    ytraj = new Tdata[d*(td->ops->max_iter+1)];
    ntraj = new Tint[td->ops->max_iter];
    xnbrI = new Tint[td->ops->max_iter];
    traj_is_tmp = true;
  } else {
    ytraj = td->out->ytraj;
    ntraj = td->out->ntraj;
    xnbrI = td->out->xnbrI;
    traj_is_tmp = false;
  }
  R2traj = new Tdata[td->ops->max_iter];



  // Loop over the probe points assigned to this particular thread
  for (probeIndex = td->probePtIndex_start; probeIndex < td->probePtIndex_end; probeIndex++) {
    isdone = false;
    iscycling = false;
    thisY = td->out->y + d*probeIndex;  // the current probe point
    thisYend = thisY+d;
    // Put the starting location on the trajectory history
    ytrajtmp = ytraj;
    for (i = 0; i < d; i++)
      ytrajtmp[i] = thisY[i];
    // We might abort when the closest landmark (or data point, when
    // applicable) changes, so keep track of it.
    closestIndexOld = -1;
    iter = 0;
    while (!isdone) {
      // if (probeIndex == 0) {
      // 	mexPrintf("iter %d, closestIndexOld = %d, closestIndex = %d\n",iter,closestIndexOld,closestIndex);
      // }
      // Do all the preparation to find neighbors efficiently
      lm_nbrs.initialize(thisY,td->x,*(td->lminfo));
      closestDataIndex = lm_nbrs.current_point()->index;
      td->out->closestDataIndex[probeIndex] = closestDataIndex + td->lminfo->index_offset;
      // Test for early termination
      if (td->ops->terminate_mode == 'x')
	// Test via the closest data point
	closestIndex = closestDataIndex;
      else if (td->ops->terminate_mode == 'l')
	// Test via the landmark assigned to the closest data point
	closestIndex = lm_nbrs.landmark_of_current_point();
      // See whether the closest landmark/point changed
      if (closestIndexOld != closestIndex && closestIndexOld >= 0)
	break; // quit on this probe point, we've reached a stopping criterion
      closestIndexOld = closestIndex;

      expand_neighborhood_msams(td->x,td->d,thisY,lm_nbrs,td->ops,x_cum,deltax_cum,sumsq_cum,x_backtrack,n,R2,td->lminfo->index_offset,(Tint *) NULL);

      // The MSAMS criterion has been met, in one form or another
      if (lm_nbrs.is_empty()) {
	// We used the entire data set.  Put y at the mean of the whole set
	for (i = 0; i < d; i++)
	  thisY[i] = x_cum[i]/n;
	td->out->n[probeIndex] = n;
	td->out->R2[probeIndex] = R2;
	isdone = true;
      } else {
	// We satisfied the MSAMS criterion without using the whole
	// data set. Update the current y position, using the
	// backtracked-value if appropriate
	if (td->ops->backtrack)
	  for (i = 0; i < d; i++)
	    thisY[i] = x_backtrack[i]/n;
	else
	  for (i = 0; i < d; i++)
	    thisY[i] = x_cum[i]/n;
	td->out->n[probeIndex] = n;
	td->out->R2[probeIndex] = R2;

	// The test for convergence by means of closer approach to a
	// different point/landmark is done above (and will be done
	// again on the next iteration).  Here we have to test for
	// cycling, i.e., revisiting the same sequence of positions.
	// Put current point on trajectory history
	ytrajtmp = ytraj+(iter+1)*d;
	for (i = 0; i < d; i++)
	  ytrajtmp[i] = thisY[i];
	ntraj[iter] = td->out->n[probeIndex];
	R2traj[iter] = td->out->R2[probeIndex];
	xnbrI[iter] = closestDataIndex;
	// Check previous points in reverse order since it's likely
	// the cycling starts at the end of the trajectory
	for (i = iter-1, ytrajtmp = ytraj+(i+1)*d; i >= 0; i--, ytrajtmp-=d)
	  if (sqrdist(thisY,thisYend,ytrajtmp) < tol*R2/n) {
	    isdone = true;
	    iscycling = true;
	    break;
	  }
	if (iscycling && i < iter-1) {
	  // It's cycling but not stable. Find the place where the
	  // maximum density was achieved, which is indicated by the
	  // largest n.  (In case of a tie, use the smallest R2 to
	  // break the tie.)
	  int imax,istart;
	  Tint nmax;
	  Tdata R2min;
	  // By default keep the last visit, since that's the one that
	  // is least likely to be contaminated by roundoff
	  nmax = ntraj[iter];
	  imax = iter;
	  R2min = R2traj[iter];
	  closestDataIndex = xnbrI[iter];
	  istart = i;
	  // Note: don't choose the first one, that's most likely to
	  // be contaminated by roundoff errors
	  for (i++; i < iter; i++)
	    if ((ntraj[i] > nmax) || (ntraj[i] == nmax && R2traj[i] < R2min)) {
	      nmax = ntraj[i];
	      imax = i;
	      R2min = R2traj[i];
	      closestDataIndex = xnbrI[i];
	    }
	  // Copy that location to the current point. 
	  /*
	  // We actually copy the next location, as that one is the
	  // one that is guaranteed to not have roundoff errors (all
	  // points that end up at this location should include the
	  // same points, in the same order)

	  // Update: this seems like a bad idea.  Better to be just a
	  // little sloppier about the tolerance.  The problem is you
	  // might get trapped in some weird local constellation of
	  // points---better to stick with the "Olympian view".  So
	  // this bit has been commented out.
	  imax++;
	  if (imax > iter)
	    imax = iter;
	  */
	  ytrajtmp = ytraj + imax*d;
	  for (i = 0; i < d; i++)
	    thisY[i] = ytrajtmp[i];
	  td->out->n[probeIndex] = nmax;
	  td->out->R2[probeIndex] = R2min;
	  td->out->closestDataIndex[probeIndex] = closestDataIndex + td->lminfo->index_offset;
	}  // if (iscycling)
      }  // end of processing after satisfying MSAMS criterion
      iter++;
      isdone = isdone || (iter >= td->ops->max_iter);
    } // Convergence of probe point   while(!isdone)
    // Do the last bit of reporting
    td->out->n_iter[probeIndex] = iter;
    if (iter < td->ops->max_iter)
      td->out->convergedFlag[probeIndex] = 1+iscycling;
    else
      td->out->convergedFlag[probeIndex] = 0;
  } // Loop over probe points


  // Clean up
  if (traj_is_tmp) {
    delete[] ytraj;
    delete[] ntraj;
    delete[] xnbrI;
  }
  delete[] R2traj;
  return;
}

template <class Tdata,class Tint>
void expand_neighborhood_msams(const Tdata *x,int d,const Tdata *thisY,landmarked_neighbors<Tdata,Tint> &lm_nbrs,const optionStruct *ops,vector<Tdata> &x_cum,vector<Tdata> &deltax_cum,vector<Tdata> &sumsq_cum,vector<Tdata> &x_backtrack,int &n,Tdata &R2,int offset,Tint *pointIndex)
{
  // A lot of intermediate storage is passed so that we don't have to
  // allocate fresh storage on each call

  int i;
  int n_backtrack; // the number of points at which ratio first exceeded 1
  bool msams_criterion_met, msams_backtrack_met, this_backtrack_met;
  Tdata tmpR,tmpR2,R2_backtrack;
  int thisPointIndex;
  const Tdata *thisPoint;

  const Tdata *thisYend = thisY + d;
  Tdata factor2 = ops->factor;
  factor2 = factor2*factor2;
  Tdata tol = ops->convergence_thresh;

  // Initialize the cumulative variables that we'll use to check
  // the MSAMS criterion
  for (i = 0; i < d; i++) {
    x_cum[i] = 0;
    deltax_cum[i] = 0;
    sumsq_cum[i] = 0;
  }

  // Process the data points needed to satisfy the MSAMS criterion
  msams_criterion_met = false;
  msams_backtrack_met = false;
  n = 0;   // The number of points we've checked so far
  n_backtrack = 0;
  while (!msams_criterion_met && !lm_nbrs.is_empty()) {
    thisPointIndex = lm_nbrs.current_point()->index;
    if (pointIndex != NULL)
      pointIndex[n] = thisPointIndex + offset;
    n++;  // Keep track of the # of points we've tested
    // Test the MSAMS criterion at the current point (i.e., radius)
    R2 = lm_nbrs.current_point()->R2;
    thisPoint = x+d*thisPointIndex;
    this_backtrack_met = false;
    typename vector<Tdata>::iterator xIterator,dxIterator,sqIterator;
    const Tdata *dataIterator;
    if (ops->any_coordinate) {
      for (dataIterator = thisY,xIterator = x_cum.begin(),dxIterator = deltax_cum.begin(),sqIterator = sumsq_cum.begin(); dataIterator < thisYend; dataIterator++,xIterator++,dxIterator++,sqIterator++,thisPoint++) {
	// The any_coordinate == true case
	*xIterator += *thisPoint;
	tmpR = *thisPoint - *dataIterator;
	*dxIterator += tmpR;
	*sqIterator += tmpR*tmpR;
	tmpR2 = (*dxIterator)*(*dxIterator);
	msams_criterion_met = msams_criterion_met || (n >= ops->n_min && tmpR2 > factor2 * *sqIterator);
	if (tmpR2 > *sqIterator)
	  this_backtrack_met = true;
      }
    } else {
      tmpR2 = 0;
      Tdata tmpssq = 0;
      for (dataIterator = thisY,xIterator = x_cum.begin(),dxIterator = deltax_cum.begin(),sqIterator = sumsq_cum.begin(); dataIterator < thisYend; dataIterator++,xIterator++,dxIterator++,sqIterator++,thisPoint++) {
	// The any_coordinate == false case
	*xIterator += *thisPoint;
	tmpR = *thisPoint - *dataIterator;
	*dxIterator += tmpR;
	*sqIterator += tmpR*tmpR;
	tmpR2 += (*dxIterator)*(*dxIterator);
	tmpssq += *sqIterator;
      }
      msams_criterion_met = msams_criterion_met || (n >= ops->n_min && tmpR2 > factor2 * tmpssq);
      if (tmpR2 > tmpssq)
	this_backtrack_met = true;
    }
    if (this_backtrack_met) {
      if (!msams_backtrack_met) {
	// This is the first time (at least in a while) that
	// we've satisfied the backtrack criterion, so store the
	// cumulative deltax
	msams_backtrack_met = true;
	for (i = 0; i < d; i++)
	  x_backtrack[i] = x_cum[i];
	n_backtrack = n;
	R2_backtrack = R2;
      }
    } else
      msams_backtrack_met = false; //clear if pt fails this_backtrack_met
    // Advance to the next point
    lm_nbrs.next_point();
  }
  if (ops->backtrack) {
    n = n_backtrack;
    R2 = R2_backtrack;
  }
}
