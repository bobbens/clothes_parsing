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

#include "mex.h"
#include "matrix.h"
#include "assert.h"
#include <string.h>

#include "../libRegionBPP/libRegionBPP.h"
#include "../libDH/libDHStandard.h"
#include "../libSPRegion/libSPRegion.h"

//#define WITH_PARALLELINF		//if defined: parallel inference rather than parallel in the samples
//#define WITH_NOPARALLEL		//if defined: turns off any parallelism

#ifdef WITH_LATENT
#include "../libLatentSPRegion/libLatentSPRegion.h"
#endif

class mstream : public std::streambuf {
public:
protected:
	virtual std::streamsize xsputn(const char *s, std::streamsize n) {
		mexPrintf("%.*s",n,s);
		mexEvalString("drawnow;");
		return n;
	} 
	virtual int overflow(int c = EOF) {
		if (c != EOF) {
			mexPrintf("%.1s",&c);
			mexEvalString("drawnow;");
		}
		return 1;
	}
};

void FillLHS(int nlhs, mxArray *plhs[], int retVal) {
	for(int k=0;k<nlhs;++k) {
		plhs[k] = mxCreateDoubleMatrix(1, 1, mxREAL);
		double* output = mxGetPr(plhs[k]);
		*output = double(retVal);
	}
}

template <class T>
class ComputationContainer {
private:
#ifdef WITH_PARALLELINF
#ifndef WITH_NOPARALLEL
	typedef libRegionBPP<T,true> SOLVER_TYPE;
#else
	typedef libRegionBPP<T,false> SOLVER_TYPE;
#endif
	typedef libDHStandard<T,false,SOLVER_TYPE> DH_TYPE;
#else
	typedef libRegionBPP<T,false> SOLVER_TYPE;
#ifndef WITH_NOPARALLEL
	typedef libDHStandard<T,true,SOLVER_TYPE> DH_TYPE;
#else
	typedef libDHStandard<T,false,SOLVER_TYPE> DH_TYPE;
#endif
#endif

	typedef libSPRegion<T,libDHBase<T> > SP_TYPE;
public:
	ComputationContainer() {}
	~ComputationContainer() {}

	template <class B>
	void extractField(const mxArray *pm, const char* name, B* defVal) {
		const mxArray *val = mxGetField(pm, 0, name);
		if(val!=NULL) {
			if(mxGetClassID(val)==mxDOUBLE_CLASS) {
				*defVal = B(*((double*)mxGetData(val)));
			} else if(mxGetClassID(val)==mxSINGLE_CLASS) {
				*defVal = B(*((float*)mxGetData(val)));
			}
		}
	}

	int ProcessSample(const mxArray* sample_x, typename DH_TYPE::GRAPH_TYPE* graph, std::vector<size_t>* obs, int LearnORPredict) {
		const mxArray *val = mxGetField(sample_x, 0, "VariableCardinalities");
		if(val==NULL) {
			std::cout << "VariableCardinalities missing." << std::endl;
			return -1;
		}

		size_t num_total_vars = mxGetNumberOfElements(val);
		if( (mxGetClassID(val)==mxDOUBLE_CLASS && sizeof(T)==4) || (mxGetClassID(val)==mxSINGLE_CLASS && sizeof(T)==8) ) {
			std::cout << "Wrong data type of variable cardinalities." << std::endl;
			return -1;
		}
		T* VarSizes = (T*)mxGetData(val);

		val = mxGetField(sample_x, 0, "Regions");
		if(val==NULL) {
			std::cout << "Regions missing." << std::endl;
			return -1;
		}
		size_t numRegions = mxGetNumberOfElements(val);

		std::vector<std::vector<size_t> > parentsIDs;
		std::map<size_t,size_t> FileIDToGraphID;

		for(size_t k=0;k<numRegions;++k) {
			typename SOLVER_TYPE::MPNode* r = new typename SOLVER_TYPE::MPNode;
			FileIDToGraphID[k] = graph->size();

			const mxArray* region = mxGetCell(val, k);

			r->c_r = *((T*)mxGetData(mxGetField(region, 0, "c_r")));
			r->flag = int(k);
			r->SR = NULL;

			const mxArray* regionVars = mxGetField(region, 0, "VariableIndices");
			if(regionVars==NULL) {
				std::cout << "VariableIndices missing." << std::endl;
				return -1;
			}
			r->num_variables = mxGetNumberOfElements(regionVars);

			T* regionVarPtr = (T*)mxGetData(regionVars);
			r->var_ix = new size_t[r->num_variables];
			r->num_states = new size_t[r->num_variables];
			r->cum_num_states = new size_t[r->num_variables+1];
			r->cum_num_states[0] = 1;

			for(size_t v=0;v<r->num_variables;++v) {
				r->var_ix[v] = size_t(regionVarPtr[v])-1;
				r->num_states[v] = size_t(VarSizes[r->var_ix[v]]);
				r->cum_num_states[v+1] = r->cum_num_states[v]*r->num_states[v];
			}

			size_t potSize = r->cum_num_states[r->num_variables];
			r->pot = new typename SOLVER_TYPE::Data;

			const mxArray* regionFeatures = mxGetField(region, 0, "Features");
			if(regionFeatures==NULL) {
				//check whether r and pot are given directly
				const mxArray* kIndex = mxGetField(region, 0, "r");
				const mxArray* potMatrix = mxGetField(region, 0, "pot");
				if(kIndex==NULL || potMatrix==NULL) {
					std::cout << "Features missing." << std::endl;
					return -1;
				} else {
					size_t numFeatures = mxGetNumberOfElements(kIndex);
					size_t numFeatures1 = mxGetN(potMatrix);
					if(numFeatures!=numFeatures1) {
						std::cout << "Dimensions of r and pot do not match." << std::endl;
						return -1;
					}
					r->featureIDs.assign(numFeatures, 0);
					r->potOrig.assign(numFeatures, NULL);
					T* data = (T*)mxGetData(potMatrix);
					for(size_t rx=0;rx<numFeatures;++rx) {
						r->featureIDs[rx] = size_t(*((T*)mxGetData(kIndex) + rx))-1;
						r->potOrig[rx] = new typename SOLVER_TYPE::Data;

						if(!mxIsSparse(potMatrix)) {
							r->potOrig[rx]->ResetDataNoDataModify(data+rx*potSize,potSize);
						} else {
							mwIndex* jcPtr = mxGetJc(potMatrix);
							mwIndex* irPtr = mxGetIr(potMatrix);							
							size_t Nrows = jcPtr[rx+1] - jcPtr[rx];
							r->potOrig[rx]->ResetDataNew(data+jcPtr[rx],irPtr+jcPtr[rx],Nrows,false);
						}
					}
				}
			} else {
				size_t numFeatures = mxGetNumberOfElements(regionFeatures);
				r->featureIDs.assign(numFeatures, 0);
				r->potOrig.assign(numFeatures, NULL);
				for(size_t rx=0;rx<numFeatures;++rx) {
					const mxArray* feature = mxGetCell(regionFeatures, rx);
					r->featureIDs[rx] = size_t(*((T*)mxGetData(mxGetField(feature, 0, "r"))))-1;
					r->potOrig[rx] = new typename SOLVER_TYPE::Data;
					const mxArray* mxPotArr = mxGetField(feature, 0, "pot");
					if(mxPotArr==NULL) {
						std::cout << "pot missing." << std::endl;
						return -1;
					}
					if(!mxIsSparse(mxPotArr)) {
						r->potOrig[rx]->ResetDataNoDataModify((T*)mxGetData(mxPotArr),potSize);
					} else {
						mwIndex* jcPtr = mxGetJc(mxPotArr);
						mwIndex* irPtr = mxGetIr(mxPotArr);
						T* data = (T*)mxGetData(mxPotArr);
						size_t NColumns = mxGetN(mxPotArr);
						if(NColumns!=1) {
							std::cout << "Number of columns of sparse matrix should be 1 while it is " << NColumns << "." << std::endl;
							return -1;
						}
						size_t Nrows = jcPtr[1] - jcPtr[0];
						r->potOrig[rx]->ResetDataNew(data,irPtr,Nrows,false);
					}
				}
			}
			
			regionFeatures = mxGetField(region, 0, "Loss");
			if(regionFeatures==NULL) {
				std::cout << "Loss missing." << std::endl;
				return -1;
			}
			r->loss = new typename SOLVER_TYPE::Data;
			if(!mxIsSparse(regionFeatures)) {
				r->loss->ResetDataNoDataModify((T*)mxGetData(regionFeatures),potSize);
			} else {
				mwIndex* jcPtr = mxGetJc(regionFeatures);
				mwIndex* irPtr = mxGetIr(regionFeatures);
				T* data = (T*)mxGetData(regionFeatures);
				size_t NColumns = mxGetN(regionFeatures);
				if(NColumns!=1) {
					std::cout << "Number of columns of sparse matrix should be 1 while it is " << NColumns << "." << std::endl;
					return -1;
				}
				size_t Nrows = jcPtr[1] - jcPtr[0];
				r->loss->ResetDataNew(data,irPtr,Nrows,false);
			}

			const mxArray* parents = mxGetField(region, 0, "Parents");
			if(parents==NULL) {
				std::cout << "Parents missing." << std::endl;
				return -1;
			}
			T* parentIXPtr = (T*)mxGetData(parents);
			size_t numParents = mxGetNumberOfElements(parents);
			std::vector<size_t> ParentConnections(numParents, 0);
			for(size_t p=0;p!=numParents;++p) {
				ParentConnections[p] = size_t(parentIXPtr[p])-1;
			}
			parentsIDs.push_back(ParentConnections);

			graph->push_back(r);
		}

		const mxArray *MXObs = mxGetField(sample_x, 0, "Observation");
		if(MXObs==NULL && LearnORPredict<2) {
			std::cout << "Observation missing." << std::endl;
			return -1;
		}
		if(MXObs!=NULL) {
			T* MXObsPtr = (T*)mxGetData(MXObs);
			obs->assign(mxGetNumberOfElements(MXObs), 0);
			for(size_t v=0;v<obs->size();++v) {
				obs->at(v) = size_t(MXObsPtr[v])-1;
			}
		} else {
			obs->assign(num_total_vars,0);
		}

		for(size_t k=0,k_e=parentsIDs.size();k!=k_e;++k) {
			typename SOLVER_TYPE::MPNode* child_node = graph->at(k);
			size_t potSize = child_node->cum_num_states[child_node->num_variables];

			for(size_t m=0,m_e=parentsIDs[k].size();m!=m_e;++m) {
				typename SOLVER_TYPE::MPNode* parent_node = graph->at(FileIDToGraphID[parentsIDs[k][m]]);

				typename SOLVER_TYPE::MsgContainer ptmp;
				ptmp.node = parent_node;
				ptmp.lambda = new T[potSize];
				memset(ptmp.lambda, 0, sizeof(T)*potSize);
				child_node->Parents.push_back(ptmp);

				typename SOLVER_TYPE::MsgContainer ctmp;
				ctmp.node = child_node;
				ctmp.lambda = ptmp.lambda;
				parent_node->Childs.push_back(ctmp);
			}
		}

		return 0;
	}

	int ClearGraph(typename std::vector<typename DH_TYPE::GRAPH_TYPE>* graph) {
		std::set<typename SOLVER_TYPE::Data*> ToBeDeleted;
		for(typename std::vector<typename DH_TYPE::GRAPH_TYPE>::iterator x=graph->begin(),x_e=graph->end();x!=x_e;++x) {
			for(typename std::vector<typename SOLVER_TYPE::MPNode*>::iterator rit=x->begin(),rit_e=x->end();rit!=rit_e;++rit) {
				typename SOLVER_TYPE::MPNode* r = *rit;
				if(r->cum_num_states!=NULL) {
					delete [] r->cum_num_states;
				}
				if(r->var_ix!=NULL) {
					delete [] r->var_ix;
				}
				if(r->num_states!=NULL) {
					delete [] r->num_states;
				}
				if(r->pot!=NULL) {
					ToBeDeleted.insert(r->pot);
				}
				r->featureIDs.clear();
				for(typename std::vector<typename SOLVER_TYPE::Data*>::iterator iter=r->potOrig.begin(),iter_e=r->potOrig.end();iter!=iter_e;++iter) {
					delete *iter;
				}
				r->potOrig.clear();
				if(r->loss!=NULL) {
					delete r->loss;
				}
				for(typename std::vector<typename SOLVER_TYPE::MsgContainer>::iterator ci=r->Childs.begin(),ci_e=r->Childs.end();ci!=ci_e;++ci) {
					if(ci->lambda!=NULL) {
						delete [] ci->lambda;
					}
				}
				delete r;
			}
		}
		return 0;
	}

	int OutputSettings(struct DH_TYPE::SPParameters& params) {
		std::cout << "Settings:" << std::endl;
		std::cout << "  ArmijoIterations:     " << params.ArmijoIterations << std::endl;
		std::cout << "  BetheCountingNumbers: " << params.BetheCountingNumbers << std::endl;
		std::cout << "  C:                    " << params.C << std::endl;
		std::cout << "  CCCPIterations:       " << params.CCCPIterations << std::endl;
		std::cout << "  CRFEraseMessages:     " << params.CRFEraseMessages << std::endl;
		std::cout << "  CRFIterations:        " << params.CRFIterations << std::endl;
		std::cout << "  CRFMPIterations:      " << params.CRFMPIterations << std::endl;
		//std::cout << "  CRFOuterExchange:     " << params.CRFOuterExchange << std::endl;
		//std::cout << "  CRFOuterIterations:   " << params.CRFOuterIterations << std::endl;
		std::cout << "  epsilon:              " << params.epsilon << std::endl;
		std::cout << "  MPIterations:         " << params.MPIterations << std::endl;
		std::cout << "  p:                    " << params.p << std::endl;
		std::cout << "  ReuseMessagesForF:    " << params.ReuseMessagesForF << std::endl;
		std::cout << "  Verbosity:            " << params.Verbosity << std::endl;
		return 0;
	}

	void PerformOperation(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], int accuracy, int LearnORPredict) {
		struct DH_TYPE::SPParameters SPSettings;
		if(nrhs>1) {
			extractField(prhs[1], "ArmijoIterations", &SPSettings.ArmijoIterations);
			extractField(prhs[1], "BetheCountingNumbers", &SPSettings.BetheCountingNumbers);
			extractField(prhs[1], "C" , &SPSettings.C);
			extractField(prhs[1], "CCCPIterations", &SPSettings.CCCPIterations);
			extractField(prhs[1], "CRFEraseMessages", &SPSettings.CRFEraseMessages);
			extractField(prhs[1], "CRFIterations", &SPSettings.CRFIterations);
			extractField(prhs[1], "CRFMPIterations", &SPSettings.CRFMPIterations);
			extractField(prhs[1], "epsilon", &SPSettings.epsilon);
			extractField(prhs[1], "MPIterations", &SPSettings.MPIterations);
			extractField(prhs[1], "p", &SPSettings.p);
			extractField(prhs[1], "ReuseMessagesForF", &SPSettings.ReuseMessagesForF);
			extractField(prhs[1], "Verbosity", &SPSettings.Verbosity);
		}

		OutputSettings(SPSettings);

		size_t numSamples = mxGetNumberOfElements(prhs[0]);
		std::vector<typename DH_TYPE::GRAPH_TYPE> graph(numSamples);
		std::vector<std::vector<size_t> > obs(numSamples);
		for(size_t x=0;x<numSamples;++x) {
			const mxArray *sample_x = mxGetCell(prhs[0], x);
			if(sample_x==NULL) {
				continue;
			}

			int retVal =  ProcessSample(sample_x, &graph[x], &obs[x], LearnORPredict);

			if(retVal!=0) {
				FillLHS(nlhs,plhs,retVal);
				ClearGraph(&graph);
				return;
			}
		}

#ifdef WITH_LATENT
		LatentStructuredPrediction<T,SP_TYPE,DH_TYPE> SP(&SPSettings, &graph, &obs);
#else
		StructuredPrediction<T,SP_TYPE,DH_TYPE> SP(&SPSettings, &graph, &obs);
#endif

		std::vector<T> w(SP.GetNumFeatures(), T(1.0));

		if(nrhs>3) {
			size_t numWeights = mxGetNumberOfElements(prhs[3]);
			if(numWeights==w.size()) {
				memcpy((char*)&w[0], (char*)mxGetData(prhs[3]), numWeights*sizeof(T));
			}
		}

		std::vector<std::map<size_t,std::vector<T> > > Beliefs;

		if(LearnORPredict<2) {
			SP.Iterate(&w);		
		}

		if(LearnORPredict>0) {
			SP.Predict(&w, Beliefs);	
		}

		if(nlhs>0) {
			mwSize tmpbuffer[1];
			tmpbuffer[0] = w.size();
			if(accuracy==4) {
				plhs[0] = mxCreateNumericArray(1, (mwSize*)&tmpbuffer, mxSINGLE_CLASS, mxREAL);
			} else {
				plhs[0] = mxCreateNumericArray(1, (mwSize*)&tmpbuffer, mxDOUBLE_CLASS, mxREAL);
			}
			T* output = (T*)mxGetData(plhs[0]);
			memcpy(output, &w[0], w.size()*sizeof(T));
		}

		if(nlhs>1) {
			mwSize dims = mwSize(Beliefs.size());
			const char *fName[] = {"Regions"};
			const char *fName2[] = {"ID","Beliefs"};
			plhs[1] = mxCreateCellArray(1, &dims);

			size_t cnt = 0;
			for(typename std::vector<std::map<size_t,std::vector<T> > >::iterator x=Beliefs.begin(),x_e=Beliefs.end();x!=x_e;++x, ++cnt) {
				mxArray* sample = mxCreateStructMatrix(1,1,1,fName);
				mwSize numRegions = mwSize(x->size());
				mxArray* beliefs = mxCreateCellArray(1, &numRegions);

				mwIndex ix = 0;
				for(typename std::map<size_t, std::vector<T> >::iterator bel=x->begin(),bel_e=x->end();bel!=bel_e;++bel, ++ix) {
					mxArray* bel_r = mxCreateStructMatrix(1,1,2,fName2);

					mwSize tmpbuffer[1];
					tmpbuffer[0] = mwSize(bel->second.size());
					mxArray* mxbel;
					if(accuracy==4) {
						mxbel = mxCreateNumericArray(1,(mwSize*)&tmpbuffer,mxSINGLE_CLASS,mxREAL);
					} else {
						mxbel = mxCreateNumericArray(1,(mwSize*)&tmpbuffer,mxDOUBLE_CLASS,mxREAL);
					}
					T* output = (T*)mxGetData(mxbel);
					memcpy(output, &bel->second[0], bel->second.size()*sizeof(T));

					mxArray* ID = mxCreateDoubleMatrix(1,1,mxREAL);
					double* outputID = mxGetPr(ID);
					*outputID = double(bel->first);

					mxSetField(bel_r,0,"Beliefs",mxbel);
					mxSetField(bel_r,0,"ID",ID);

					mxSetCell(beliefs, ix, bel_r);
				}

				mxSetField(sample, 0, "Regions", beliefs);
				mxSetCell(plhs[1], cnt, sample);
			}
		}

		ClearGraph(&graph);
	}
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	mstream mout;
	std::streambuf *outbuf = std::cout.rdbuf(&mout);

	int LearnORPredict = 0;
	int accuracy = 8;//double
	if(nrhs>2) {
		if(mxGetClassID(prhs[2])==mxDOUBLE_CLASS) {
			accuracy = 8;
			LearnORPredict = int(*((double*)mxGetData(prhs[2])));
		} else if(mxGetClassID(prhs[2])==mxSINGLE_CLASS) {
			accuracy = 4;
			LearnORPredict = int(*((float*)mxGetData(prhs[2])));
		} else {
			std::cout << "Input class type not recognized." << std::endl;
			FillLHS(nlhs, plhs, -1);
			return;
		}
	}

	if(accuracy==4) {
		ComputationContainer<float> op;
		op.PerformOperation(nlhs, plhs, nrhs, prhs, accuracy, LearnORPredict);
	} else {
		ComputationContainer<double> op;
		op.PerformOperation(nlhs, plhs, nrhs, prhs, accuracy, LearnORPredict);
	}

	std::cout.rdbuf(outbuf); 
}