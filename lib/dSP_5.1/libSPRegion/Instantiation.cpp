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
#include "../libRegionBPP/libRegionBPP.h"
#include "../libDH/libDHStandard.h"
#include "../libDHMultM/libDHSampPMultM.h"
#include "../libDHMultM/libDHModlPMultM.h"
#include "libSPRegion.h"
#include "libSPRegionMultM.h"

template class StructuredPrediction<double,libSPRegion<double,libDHBase<double> >,libDHStandard<double,true,libRegionBPP<double,false> > >;
template class StructuredPrediction<float,libSPRegion<float,libDHBase<float> >,libDHStandard<float,true,libRegionBPP<float,false> > >;
template class StructuredPrediction<double,libSPRegion<double,libDHBase<double> >,libDHStandard<double,false,libRegionBPP<double,false> > >;
template class StructuredPrediction<float,libSPRegion<float,libDHBase<float> >,libDHStandard<float,false,libRegionBPP<float,false> > >;
template class StructuredPrediction<double,libSPRegion<double,libDHBase<double> >,libDHStandard<double,false,libRegionBPP<double,true> > >;
template class StructuredPrediction<float,libSPRegion<float,libDHBase<float> >,libDHStandard<float,false,libRegionBPP<float,true> > >;

template class StructuredPrediction<double,libSPRegionMultM<double,libDHBaseMultM<double> >,libDHSampPMultM<double,libRegionBPP<double,false> > >;
template class StructuredPrediction<float,libSPRegionMultM<float,libDHBaseMultM<float> >,libDHSampPMultM<float,libRegionBPP<float,false> > >;
template class StructuredPrediction<double,libSPRegionMultM<double,libDHBaseMultM<double> >,libDHModlPMultM<double,libRegionBPP<double,true> > >;
template class StructuredPrediction<float,libSPRegionMultM<float,libDHBaseMultM<float> >,libDHModlPMultM<float,libRegionBPP<float,true> > >;