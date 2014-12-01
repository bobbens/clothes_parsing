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
#ifndef __LIB__DHSTANDARD__
#define __LIB__DHSTANDARD__

#include <vector>
#include <string>

#include "libDHBase.h"
#include "../libRegionBPP/libRegionBPP.h"

template <class T, bool B, class S>
class libDHStandard : public libDHBase<T> {
	template <class Y> struct Sample {
		std::vector<Y>* x;
	};
	template<bool U> struct i2t {};
public:
	typedef typename std::vector<typename S::MPNode*> GRAPH_TYPE;
private:
	struct libDHBase<T>::SPParameters* DHParams;
	size_t NumFeatures;

	bool isMyGraph;
	bool isMyObs;
	struct Sample<GRAPH_TYPE>* phi;
	std::vector<std::vector<size_t> >* y;
	std::vector<S*> rBP;
	std::vector<T> ObservedEmpiricalMean;

	int ReadRegionFile(const char* fn, std::vector<typename S::MPNode*>& Graph, bool WithFeatureID, T (*conversion)(T) = NULL);
	int ReadRegionFileBinary(const char* fn, std::vector<typename S::MPNode*>& Graph, bool WithFeatureID, T (*conversion)(T));
	std::string ReplaceExtension(std::string orig, std::string extension);
	int ReadObservationFile(const char* fn, std::vector<size_t>& y);

	void GetAncestors(std::set<typename S::MPNode*>& Ancestors, typename S::MPNode* r);

	int ClearFeature(std::vector<typename S::MPNode*>& Graph);
	int Marginalize(std::vector<GRAPH_TYPE>* graph);
	T ComputeFandGSpecial(i2t<true>, std::vector<T> *theta, std::vector<T> *grad, int ToReturn, size_t MPIterations);
	T ComputeFandGSpecial(i2t<false>, std::vector<T> *theta, std::vector<T> *grad, int ToReturn, size_t MPIterations);
	int PredictSpecial(i2t<true>, std::vector<T> *theta, std::vector<std::map<size_t,std::vector<T> > >& Beliefs);
	int PredictSpecial(i2t<false>, std::vector<T> *theta, std::vector<std::map<size_t,std::vector<T> > >& Beliefs);
public:
	int DHClear();
	int ComputeEmpiricalMean(std::vector<T>* empMean);
	T ComputeEmpiricalMean(std::vector<T>* empMean, std::vector<T>* theta);
	int ClearMessages();
	T ComputeFandG(std::vector<T> *theta, std::vector<T> *empiricalMeans, std::vector<T> *grad, int flag);
public:
	libDHStandard(struct libDHBase<T>::SPParameters* params);
	libDHStandard(struct libDHBase<T>::SPParameters* params, std::vector<GRAPH_TYPE>* graph, std::vector<std::vector<size_t> >* obs, bool marginalize = false);
	~libDHStandard();

	size_t GetNumFeatures();
	int Predict(std::vector<T> *theta, std::vector<std::map<size_t,std::vector<T> > >& Beliefs);
	std::vector<typename libDHStandard<T,B,S>::GRAPH_TYPE>* GetGraph();
	std::vector<std::vector<size_t> >* GetObservation();
};

#endif