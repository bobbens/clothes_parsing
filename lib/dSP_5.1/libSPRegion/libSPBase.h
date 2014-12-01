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
#ifndef __LIB__SPBASE__
#define __LIB__SPBASE__

#include <cstddef>
#include <vector>

template <class T, class DHB>
class libSPBase {
protected:
	virtual int SPClear() = 0;
	DHB* DataHandler;
public:
	virtual ~libSPBase() {};
	virtual T Iterate(std::vector<T>* w, std::vector<T>* empMean) = 0;
};

template <class T, class SP, class DH>
class StructuredPrediction : public SP, public DH {
public:
	StructuredPrediction(struct DH::SPParameters* spp) : SP(spp), DH(spp) {SP::DataHandler = this;}
	StructuredPrediction(struct DH::SPParameters* spp, std::vector<typename DH::GRAPH_TYPE>* graph, std::vector<std::vector<size_t> >* obs) : SP(spp), DH(spp, graph, obs) {SP::DataHandler = this;}
	
	int Clear() {
		DH::DHClear();
		SP::SPClear();
		return 0;
	}
};

#endif