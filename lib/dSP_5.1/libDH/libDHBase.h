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
#ifndef __LIB_DHBASE__
#define __LIB_DHBASE__

#include <vector>
#include <map>
#include <string>

template <class T, class DHB> class libSPRegion;

template <class T>
class libDHBase {
	friend class libSPRegion<T, libDHBase<T> >;
private:
	virtual int ComputeEmpiricalMean(std::vector<T>* empMean) = 0;
	virtual int ClearMessages() = 0;
	virtual T ComputeFandG(std::vector<T> *theta, std::vector<T> *empiricalMeans, std::vector<T> *grad, int flag) = 0;
public:
	virtual int DHClear() = 0;
public:
	virtual size_t GetNumFeatures() = 0;
	virtual int Predict(std::vector<T> *theta, std::vector<std::map<size_t,std::vector<T> > >& Beliefs) = 0;
	virtual ~libDHBase() {}

	struct SPParameters {
		T epsilon;
		T p;
		T C;
		size_t CRFMPIterations;
		int ReuseMessagesForF;
		size_t CRFIterations;
		size_t ArmijoIterations;
		int CRFEraseMessages;
		int BetheCountingNumbers;
		size_t CCCPIterations;
		size_t MPIterations;
		int ReadBinary;
		int Verbosity;
		std::vector<std::string> FeatureFiles;
		SPParameters() : epsilon(T(1.0)), p(T(2.0)), C(T(1.0)), CRFMPIterations(1), ReuseMessagesForF(1), CRFIterations(10), ArmijoIterations(50), CRFEraseMessages(0), BetheCountingNumbers(0), CCCPIterations(1), MPIterations(50), ReadBinary(0), Verbosity(0) {;}
	};
};

#endif