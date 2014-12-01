/*  This file is part of distributed Structured Prediction (dSP) - http://www.alexander-schwing.de/
 *
 *  distributed Structured Prediction (dSP) is free software: you can
 *  redistribute it and/or modify it under the terms of the GNU General
 *  Public License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  distributed Structured Prediction (dSP) is distributed in the hope
 *  that it will be useful, but WITHOUT ANY WARRANTY; without even the
 *  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE. See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with distributed Structured Prediction (dSP).
 *  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright (C) 2010-2013  Alexander G. Schwing  [http://www.alexander-schwing.de/]
 */

//Author: Alexander G. Schwing

/*  ON A PERSONAL NOTE: I spent a significant amount of time to go through
 *  both, theoretical justifications and coding of this framework.
 *  I hope the package is useful for your task. Any requests, 
 *  feedback, donations and support would be greatly appreciated.
 *  Thank you for contacting me!
 */
#ifndef __LIB_LATENTSPREGION__
#define __LIB_LATENTSPREGION__

#include <vector>
#include <map>

#include "../libSPRegion/libSPBase.h"

template <class T, class SP, class DH>
class LatentStructuredPrediction {
	StructuredPrediction<T,SP,DH>* StructPred;
	DH* DataHandler;
public:
	typename DH::SPParameters* SPParams;
public:
	LatentStructuredPrediction(struct DH::SPParameters* spp);
	LatentStructuredPrediction(struct DH::SPParameters* spp, std::vector<typename DH::GRAPH_TYPE>* graph, std::vector<std::vector<size_t> >* obs);
	
	T Iterate(std::vector<T>* theta);
	int Clear();
	int Predict(std::vector<T> *theta, std::vector<std::map<size_t,std::vector<T> > >& Beliefs);
	size_t GetNumFeatures();
};

#endif