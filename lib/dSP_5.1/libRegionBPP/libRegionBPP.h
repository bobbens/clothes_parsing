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
#ifndef __LIB_REGIONBPP__
#define __LIB_REGIONBPP__

#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <assert.h>

#include "DataContainers.h"
#include "../GeneralDefinitions.h"

template <class T, bool B>
class libRegionBPP {
	template<bool U> struct i2t {};
public:
	struct PotentialStruct {
		T* pot;
		size_t* num_states;
		size_t num_variables;
		PotentialStruct() : pot(NULL), num_states(NULL), num_variables(-1) {;}
		PotentialStruct(T* p, size_t* num_s, size_t num_v) : pot(p), num_states(num_s), num_variables(num_v) {;}
	};

	typedef DataContainer<T> Data;

	struct MPNode;

	struct SharedRegion {
		struct MPNode* Region;
		T* nu;
		size_t cumSize;
		size_t numMach;
	};

	struct MPRegion {
		int flag;
		T c_r;
		std::vector<size_t> featureIDs;
		size_t num_variables;
		size_t* num_states;
		size_t* cum_num_states;
		size_t* var_ix;
		std::vector<Data*> potOrig;
		Data* loss;
		Data* pot;
		T* bel;
		struct SharedRegion* SR;
		size_t Color;
		MPRegion() : SR(NULL) {};
	};

	struct MsgContainer {
		T* lambda;
		struct MPNode* node;
	};

	struct MPNode : MPRegion {
		std::vector<MsgContainer> Parents;
		std::vector<MsgContainer> Childs;
	};

	std::vector<MPNode*>* Graph;
private:
	T factorThreshold;
	T agreementThreshold;
	T epsilon;
	bool isMyGraph;
	std::vector<std::vector<MPNode*> > ColorToNodes;
	std::vector<MPNode*> ProcessingOrder;

	void ComputeLambdaRtoP(MPNode* r);
	T ComputeMuPToR(MPNode *p, MPNode *r, std::map<size_t,size_t>& x_r);

	T ComputePhiHat(MPNode* r, size_t x_r, size_t numVars);
	T ComputePrimal();
	int ComputeRegionBelief(MPNode *r, T* mem);
	T Entropy(T* mem, size_t numEL);

	int ConvertToGraph(std::map<size_t, std::pair<struct PotentialStruct, std::set< size_t > > >& LocalBPPot, std::map<size_t, std::pair<struct PotentialStruct, std::set< size_t > > >& FactorBPPot);
	int GraphColoring(i2t<true>);
	int GraphColoring(i2t<false>);
	int TreeDetection();

	int OutputMessages();
	void IterateTree();
	void Iterate(i2t<true>);
	void Iterate(i2t<false>);
	T ComputeDual(i2t<true>);
	T ComputeDual(i2t<false>);
	T ComputeDualDo(MPNode* r);
	int GetResult(i2t<true>, std::map<size_t,std::vector<T> >& Beliefs);
	int GetResult(i2t<false>, std::map<size_t,std::vector<T> >& Beliefs);
public:
	libRegionBPP(std::vector<MPNode*>* Graph, T epsilon);
	libRegionBPP(std::map<size_t, std::pair<struct PotentialStruct, std::set< size_t > > >& LocalBPPot, std::map<size_t, std::pair<struct PotentialStruct, std::set< size_t > > >& FactorBPPot, T epsilon);
	~libRegionBPP();

	void SetEpsilon(T epsilon) {this->epsilon = epsilon;}
	int Initialize();
	T RunMP(size_t MaxIter, T pdGapThr, int ToReturn = 0, int verbose = 2, CPrecisionTimer* CTmr = NULL);
	int GetResult(std::map<size_t,std::vector<T> >& Beliefs);

	T ComputeDual();
	int ComputePrimalWithAgreement(T *primal, T *primalAgree, T* EntropySum);
};

#endif