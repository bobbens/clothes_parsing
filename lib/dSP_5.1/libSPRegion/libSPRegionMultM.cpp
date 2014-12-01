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
#include <mpi.h>

#include "../libDHMultM/libDHBaseMultM.h"
#include "libSPRegionMultM.h"
#include "../GeneralDefinitions.h"

#define sign2(x) (( (x) > 0 ) - ( (x) < 0 ))

template <class T, class DHB>
libSPRegionMultM<T,DHB>::libSPRegionMultM(typename DHB::SPParameters* params) {
	SPParams = params;
	empiricalMean = NULL;
	myEmpiricalMean = false;
}

template <class T, class DHB>
libSPRegionMultM<T, DHB>::~libSPRegionMultM() {
	SPClear();
}

template <class T, class DHB>
int libSPRegionMultM<T, DHB>::SPClear() {
	if(empiricalMean!=NULL && myEmpiricalMean) {
		empiricalMean->clear();
		std::vector<T>().swap(*empiricalMean);
		delete empiricalMean;
		empiricalMean = NULL;
	}
	return 0;
}

template <class T, class DHB>
T libSPRegionMultM<T,DHB>::Iterate(std::vector<T>* w, std::vector<T>* empMean) {
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

	CPrecisionTimer* tmr = NULL;
	if(libSPBase<T,DHB>::DataHandler->GetRank()==0) {
		tmr = new CPrecisionTimer;
		tmr->Start();
	}

	T oldPrimalValue = 0;
	assert(w->size()==NumFeatures);
	std::vector<T> grad(NumFeatures);
	for(size_t iteration=0,iteration_e=SPParams->CRFIterations;iteration<iteration_e;++iteration) {
		if(libSPBase<T,DHB>::DataHandler->GetRank()==0) {
			std::cout << "I: " << iteration;
		}
		oldPrimalValue = ArmijoStep(w, empiricalMean, grad, iteration);

		if(SPParams->CRFEraseMessages) {
			libSPBase<T,DHB>::DataHandler->ClearMessages();
		}

		if(libSPBase<T,DHB>::DataHandler->GetRank()==0) {
			std::cout << " T: " << tmr->Stop() << std::endl;
		}
	}

	if(tmr!=NULL) {
		delete tmr;
	}

	return oldPrimalValue;
}

template <class T, class DHB>
T libSPRegionMultM<T,DHB>::ArmijoStep(std::vector<T> *theta, std::vector<T> *empMean, std::vector<T>& grad, size_t iteration) {
	std::vector<T> newTheta(*theta);
	grad.assign(NumFeatures, T(0.0));

	int flag = 0;
	if((iteration+1)%SPParams->CRFOuterExchange==0) {
		flag = 1;
	}

	T oldPrimalValue = libSPBase<T,DHB>::DataHandler->ComputeFandG(&newTheta, empMean, &grad, flag);

	++funEvals;

	T stepsize = T(1.0);

	for(size_t armijoIteration=0;armijoIteration<SPParams->ArmijoIterations;++armijoIteration) {
		if(libSPBase<T,DHB>::DataHandler->GetRank()==0) {
			for(size_t r = 0;r<NumFeatures;++r) {
				newTheta[r] = (*theta)[r] - stepsize*grad[r];
			}	
		}
		if(sizeof(T)==8) {
			MPI::COMM_WORLD.Bcast(&newTheta[0], int(NumFeatures), MPI::DOUBLE, 0);
		} else if(sizeof(T)==4) {
			MPI::COMM_WORLD.Bcast(&newTheta[0], int(NumFeatures), MPI::FLOAT, 0);
		}

		if(SPParams->ReuseMessagesForF==0) {
			libSPBase<T,DHB>::DataHandler->ClearMessages();
		}
		T primalValue = libSPBase<T,DHB>::DataHandler->ComputeFandG(&newTheta, empMean, NULL, 0);
		++funEvals;

		if(libSPBase<T,DHB>::DataHandler->GetRank()==0) {
			int ContinueFlag = 1;
			if(primalValue<oldPrimalValue) {
				for(size_t r=0;r<NumFeatures;++r) {
					(*theta)[r] = newTheta[r];
				}
				oldPrimalValue = primalValue;
				ContinueFlag = 0;
				MPI::COMM_WORLD.Bcast(&ContinueFlag, 1, MPI::INT, 0);
				if(sizeof(T)==8) {
					MPI::COMM_WORLD.Bcast(&(*theta)[0], int(NumFeatures), MPI::DOUBLE, 0);
				} else if(sizeof(T)==4) {
					MPI::COMM_WORLD.Bcast(&(*theta)[0], int(NumFeatures), MPI::FLOAT, 0);
				}
				break;
			} else {
				MPI::COMM_WORLD.Bcast(&ContinueFlag, 1, MPI::INT, 0);
			}
			stepsize /= 2;
		} else {
			int ContinueFlag = 1;
			MPI::COMM_WORLD.Bcast(&ContinueFlag, 1, MPI::INT, 0);
			if(ContinueFlag==0) {
				if(sizeof(T)==8) {
					MPI::COMM_WORLD.Bcast(&(*theta)[0], int(NumFeatures), MPI::DOUBLE, 0);
				} else if(sizeof(T)==4) {
					MPI::COMM_WORLD.Bcast(&(*theta)[0], int(NumFeatures), MPI::FLOAT, 0);
				}
				break;
			}
		}
	}

	if(libSPBase<T,DHB>::DataHandler->GetRank()==0) {
		std::cout << " F: " << funEvals << " S: " << stepsize << " P: " << oldPrimalValue;
	}
	return oldPrimalValue;
}

template class libSPRegionMultM<double,libDHBaseMultM<double> >;
template class libSPRegionMultM<float, libDHBaseMultM<float> >;