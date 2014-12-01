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
#include "libLatentSPRegion.h"

#include "../libSPRegion/libSPRegion.h"
#include "../libSPRegion/libSPRegionMultM.h"

#include "../libDH/libDHStandard.h"
#include "../libDHMultM/libDHSampPMultM.h"

template <class T, class SP, class DH>
LatentStructuredPrediction<T,SP,DH>::LatentStructuredPrediction(struct DH::SPParameters* spp) {
	SPParams = spp;
	StructPred = new StructuredPrediction<T,SP,DH>(spp);
	DataHandler = new DH(spp, StructPred->GetGraph(), StructPred->GetObservation(), true);
}

template <class T, class SP, class DH>
LatentStructuredPrediction<T,SP,DH>::LatentStructuredPrediction(struct DH::SPParameters* spp, std::vector<typename DH::GRAPH_TYPE>* graph, std::vector<std::vector<size_t> >* obs) {
	SPParams = spp;
	StructPred = new StructuredPrediction<T,SP,DH>(spp,graph,obs);
	DataHandler = new DH(spp, graph, obs, true);
}

template <class T, class SP, class DH>
int LatentStructuredPrediction<T,SP,DH>::Clear() {
	if(StructPred!=NULL) {
		StructPred->Clear();
		delete StructPred;
		StructPred=NULL;
	}
	if(DataHandler!=NULL) {
		DataHandler->DHClear();
		delete DataHandler;
		DataHandler = NULL;
	}
	return 0;
}

template <class T, class SP, class DH>
size_t LatentStructuredPrediction<T,SP,DH>::GetNumFeatures() {
	return DataHandler->GetNumFeatures();
}

template <class T, class SP, class DH>
T LatentStructuredPrediction<T,SP,DH>::Iterate(std::vector<T> *theta) {
	std::vector<T> empMean;
	T HCRFPrimal = T(0.0);
	for(size_t k=0;k<SPParams->CCCPIterations;++k) {
		T EntropySum = DataHandler->ComputeEmpiricalMean(&empMean, theta);
		HCRFPrimal = StructPred->Iterate(theta, &empMean) + EntropySum;
	}
	return HCRFPrimal;
}

template <class T, class SP, class DH>
int LatentStructuredPrediction<T,SP,DH>::Predict(std::vector<T> *theta, std::vector<std::map<size_t,std::vector<T> > >& Beliefs) {
	return StructPred->Predict(theta, Beliefs);
}

template class LatentStructuredPrediction<double,libSPRegion<double,libDHBase<double> >,libDHStandard<double,true,libRegionBPP<double,false> > >;
template class LatentStructuredPrediction<float,libSPRegion<float,libDHBase<float> >,libDHStandard<float,true,libRegionBPP<float,false> > >;
template class LatentStructuredPrediction<double,libSPRegion<double,libDHBase<double> >,libDHStandard<double,false,libRegionBPP<double,false> > >;
template class LatentStructuredPrediction<float,libSPRegion<float,libDHBase<float> >,libDHStandard<float,false,libRegionBPP<float,false> > >;
template class LatentStructuredPrediction<double,libSPRegion<double,libDHBase<double> >,libDHStandard<double,false,libRegionBPP<double,true> > >;
template class LatentStructuredPrediction<float,libSPRegion<float,libDHBase<float> >,libDHStandard<float,false,libRegionBPP<float,true> > >;

template class LatentStructuredPrediction<double,libSPRegionMultM<double,libDHBaseMultM<double> >,libDHSampPMultM<double,libRegionBPP<double,false> > >;
template class LatentStructuredPrediction<float,libSPRegionMultM<float,libDHBaseMultM<float> >,libDHSampPMultM<float,libRegionBPP<float,false> > >;