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
#ifndef __lib_SPREGIONMULTM__
#define __lib_SPREGIONMULTM__

#include <vector>

#include "libSPBase.h"

template<class T, class DHB>
class libSPRegionMultM : public libSPBase<T, DHB> {
public:

private:
	typename DHB::SPParameters* SPParams;

	size_t NumFeatures;
	size_t funEvals;
	std::vector<T>* empiricalMean;
	bool myEmpiricalMean;

	T ArmijoStep(std::vector<T> *theta, std::vector<T> *empMean, std::vector<T>& grad, size_t iteration);
protected:
	int SPClear();
public:
	libSPRegionMultM(typename DHB::SPParameters* params);
	~libSPRegionMultM();

	T Iterate(std::vector<T>* w, std::vector<T>* empMean = NULL);
};

#endif