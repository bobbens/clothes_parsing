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
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>
#include <string.h>
#include <assert.h>

#include <omp.h>

#include "../libDH/libDHBase.h"
#include "libSPRegion.h"
#include "../GeneralDefinitions.h"

#define sign2(x) (( (x) > 0 ) - ( (x) < 0 ))

template <class T, class DHB>
libSPRegion<T,DHB>::libSPRegion(typename DHB::SPParameters* params) {
	SPParams = params;
	empiricalMean = NULL;
	myEmpiricalMean = false;
}

template <class T, class DHB>
libSPRegion<T,DHB>::~libSPRegion() {
	SPClear();
}

template <class T, class DHB>
int libSPRegion<T,DHB>::SPClear() {
	if(empiricalMean!=NULL && myEmpiricalMean) {
		empiricalMean->clear();
		std::vector<T>().swap(*empiricalMean);
		delete empiricalMean;
		empiricalMean = NULL;
	}
	return 0;
}

template <class T, class DHB>
T libSPRegion<T, DHB>::Iterate(std::vector<T>* w, std::vector<T>* empMean) {
	NumFeatures = libSPBase<T,DHB>::DataHandler->GetNumFeatures();
	funEvals = 0;
	if(empMean==NULL) {
		myEmpiricalMean = true;
		empiricalMean = new std::vector<T>(NumFeatures, T(0.0));
		libSPBase<T,DHB>::DataHandler->ComputeEmpiricalMean(empiricalMean);
	} else {
		myEmpiricalMean = false;
		empiricalMean = empMean;
	}

	CPrecisionTimer tmr;
	tmr.Start();

	T oldPrimalValue = 0;
	assert(w->size()==NumFeatures);
	std::vector<T> grad(NumFeatures);
	for(size_t iteration=0,iteration_e=SPParams->CRFIterations;iteration<iteration_e;++iteration) {
		std::cout << "I: " << iteration;
		oldPrimalValue = ArmijoStep(w, empiricalMean, grad);

		if(SPParams->CRFEraseMessages) {
			libSPBase<T,DHB>::DataHandler->ClearMessages();
		}

		std::cout << " T: " << tmr.Stop() << std::endl;
	}

	return oldPrimalValue;
}

template <class T, class DHB>
T libSPRegion<T, DHB>::ArmijoStep(std::vector<T> *theta, std::vector<T> *empMean, std::vector<T>& grad) {
	std::vector<T> newTheta(*theta);
	grad.assign(NumFeatures, T(0.0));

	T oldPrimalValue = libSPBase<T,DHB>::DataHandler->ComputeFandG(&newTheta, empMean, &grad, 0);
	++funEvals;

	T stepsize = T(1.0);

	for(size_t armijoIteration=0;armijoIteration<SPParams->ArmijoIterations;++armijoIteration) {
		for(size_t r = 0;r<NumFeatures;++r) {
			newTheta[r] = (*theta)[r] - stepsize*grad[r];
		}

		if(SPParams->ReuseMessagesForF==0) {
			libSPBase<T,DHB>::DataHandler->ClearMessages();
		}
		T primalValue = libSPBase<T,DHB>::DataHandler->ComputeFandG(&newTheta, empMean, NULL, 0);
		++funEvals;

		if(primalValue<oldPrimalValue) {
			for(size_t r=0;r<NumFeatures;++r) {
				(*theta)[r] = newTheta[r];
			}
			oldPrimalValue = primalValue;
			break;
		}
		stepsize /= 2;
	}

	std::cout << " F: " << funEvals << " S: " << stepsize << " P: " << oldPrimalValue;
	return oldPrimalValue;
}

template class libSPRegion<double, libDHBase<double> >;
template class libSPRegion<float, libDHBase<float> >;