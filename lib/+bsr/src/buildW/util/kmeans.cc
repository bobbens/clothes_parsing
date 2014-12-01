// Copyright (C) 2002 Charless C. Fowlkes <fowlkes@eecs.berkeley.edu>
// Copyright (C) 2002 David R. Martin <dmartin@eecs.berkeley.edu>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
// 02111-1307, USA, or see http://www.gnu.org/copyleft/gpl.html.

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <assert.h>

#include "util.hh"
#include "timer.hh"
#include "random.hh"
#include "array.hh"
#include "kmeans.hh"

#define DEBUG_LEVEL 0

namespace Util
{

  // perform kmeans clustering.  multi-level implementation.
  //
  // points       - input data, numPoints rows x numDim columns
  // numPoints    - number of input data points
  // numDim       - dimensionality of input data points
  // numClusters  - number of desired clusters
  // minChange    - stop iterating if fractional change in RMS error 
  //                is less than this value
  // numIterations - maximum number of iterations
  // means        - output centroids, numClusters rows x numDim columns
  // membership   - output membership, numPoints values in [0,numClusters) 
  // multiLevel   - do multi-level?
  //
  // means and memberships must be pre-allocated.
  //
  float
  kmeans (const Array2D<float> points, const int numClusters, const float minChange, const int numIterations, 
          const bool multiLevel, Array2D<float>& means, Array1D<int>& membership)
  {
      Util::Message::startBlock("kmeans",DEBUG_LEVEL);
      const int numPoints = points.size(0);
      const int numDim = points.size(1);
      assert (numPoints >= numClusters);
      
      means.resize(numClusters,numDim);
      membership.resize(numPoints);
      membership.init(-1);

      Util::Array1D <int> count(numClusters);
      count.init (0);

      const int kmeansCoarsenRate = 3;
      const int coarseNumPoints = numPoints / kmeansCoarsenRate;
      if (!multiLevel || coarseNumPoints < numClusters)
      {
          // Not enough points to coarsen the problem.
          // Initialize membership randomly; make sure each cluster
          // has at least one point.
          Util::Array1D <int>pick (numClusters);
          Util::kOfN (numClusters, numPoints, (uint *) pick.data ());
          // Assign one random point to each cluster.
          for (int i = 0; i < numClusters; i++)
          {
            const int cluster = i;
            const int j = pick (i);
            assert (j >= 0 && j < numPoints);
            membership(j) = cluster;
            count (cluster)++;
          }
          // Assign remaining points to random clusters.
          for (int i = 0; i < numPoints; i++)
          {
            if (membership(i) != -1)
            {
                continue;
            }
            const int cluster = Util::rand->i32 (0, numClusters - 1);
            membership(i) = cluster;
            count (cluster)++;
          }
      }
      else
      {
          // Solve coarsened problem to find initial centroids.
          Util::Array1D < int >pick (coarseNumPoints);
          Util::kOfN (coarseNumPoints, numPoints, (uint *) pick.data ());
          Util::Array2D < float >coarsePoints (coarseNumPoints, numDim);
          Util::Array1D < int >coarseMembership (coarseNumPoints);
          for (int i = 0; i < coarseNumPoints; i++)
          {
            for (int k = 0; k < numDim; k++)
            {
              const int j = pick (i);
              assert (j >= 0 && j < numPoints);
              coarsePoints(i,k) = points(j,k);
            }
          }
          kmeans (coarsePoints, numClusters, minChange, numIterations,
                  multiLevel, means, coarseMembership);
          
          Util::Message::debug("assign membership from recursive call",DEBUG_LEVEL);
          for (int i = 0; i < coarseNumPoints; i++)
          {
            const int j = pick(i);
            assert (membership(j) == -1);
            const int cluster = coarseMembership(i);
            membership(j) = cluster;
            count (cluster)++;
          }
          Util::Message::debug("assign remaining points to nearest centroid",DEBUG_LEVEL);
          for (int i = 0; i < numPoints; i++)
          {
            if (membership(i) != -1)
            {
              continue;
            }
            float mindist = FLT_MAX;
            int cluster = -1;
            for (int j = 0; j < numClusters; j++)
            {
              float dist = 0;
              for (int k = 0; k < numDim; k++)
              {
                const float dif = points(i,k) - means(j,k);
                dist += dif * dif;
              }
              if (dist < mindist)
              {
                mindist = dist;
                cluster = j;
              }
            }
            //assert (cluster != -1);
            if(cluster == -1)
            {
              Util::Message::error("no nearest cluster????");   
            }
            membership(i) = cluster;
            count (cluster)++;
          }
      }

      // Make sure every cluster has at least one point.
      for (int j = 0; j < numClusters; j++)
      {
        assert (count (j) > 0);
      }

      // Compute the mean of each cluster.
      means.init(0);
      for (int i = 0; i < numPoints; i++)
      {
        const int j = membership(i);
        for (int k = 0; k < numDim; k++)
        {
          means(j,k) += points(i,k);
        }
      }
      for (int j = 0; j < numClusters; j++)
      {
        assert (count (j) > 0);
        for (int k = 0; k < numDim; k++)
        {
            means(j,k) /= (float)count(j);
        }
      }

      // iterate
      Util::Message::debug(String("kmeans: n=%d d=%d k=%d it=%d",
                                  numPoints,numDim,numClusters,numIterations),DEBUG_LEVEL);
      Util::Array1D <int> newmember (numPoints);
      float deltaerr = FLT_MAX;
      float err = FLT_MAX;
      for (int itr = 0; itr < numIterations && deltaerr > minChange; itr++)
      {
        //update stopping conditions
        float olderr = err;
        err = 0;

        // compute cluster membership
        for (int i = 0; i < numPoints; i++)
        {
          float mindist = FLT_MAX;
          int cluster = -1;
          for (int j = 0; j < numClusters; j++)
          {
            float dist = 0;
            for (int k = 0; k < numDim; k++)
            {
              const float dif = points(i,k) - means(j,k);
              dist += dif * dif;
            }
            if (dist > FLT_MAX)
              Util::Message::error(String("dist = %f",dist));

            if (dist < mindist)
            {
                mindist = dist;
                cluster = j;
            }
          }
            if(cluster == -1)
            {
              Util::Message::error("2 no nearest cluster????");   
            }
          newmember (i) = cluster;
          err += mindist;
        }
        err = sqrt (err / numPoints);

        deltaerr = (err == 0) ? 0 : ((olderr - err) / err);

        //adjust means
        for (int i = 0; i < numPoints; i++)
        {
          //move the point if necessary and adjust the means but
          //dont move if it's the last member of a cluster
          const int om = membership(i);
          const int nm = newmember(i);
          assert (om >= 0 && om < numClusters);
          assert (nm >= 0 && nm < numClusters);
          assert (count (om) > 0);
          assert (count (nm) > 0);
          if (nm != om && count (om) > 1)
          {
            const float oc = (float) count (om);
            const float nc = (float) count (nm);
            for (int k = 0; k < numDim; k++)
            {
                const float dat = points(i,k);
                means(om,k) = (oc*means(om,k) - dat) / (oc - 1);
                means(nm,k) = (nc*means(nm,k) + dat) / (nc + 1);
            }
            membership(i) = nm;
            count (om)--;
            count (nm)++;
          }
        }

      }
      Util::Message::endBlock(DEBUG_LEVEL);
      Util::Message::debug(String("RMS error = %f",err),DEBUG_LEVEL);
      return err;
  }

} //namespace Util
