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
#include <cmath>
#include <limits>

#include <string.h>
#include <omp.h>
#include <mpi.h>

#include "libDHSampPMultM.h"

#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

#define sign2(x) (( (x) > 0 ) - ( (x) < 0 ))

template <class T, class S>
libDHSampPMultM<T, S>::libDHSampPMultM(struct libDHBaseMultM<T>::SPParameters* params) {
	DHParams = params;

	MPI::Init(DHParams->argc, DHParams->argv);
	ClusterSize = MPI::COMM_WORLD.Get_size();
	Rank = MPI::COMM_WORLD.Get_rank();

	NumFeatures = 0;
	isMyGraph = true;
	isMyObs = true;
	phi = new struct Sample<GRAPH_TYPE>;
	phi->x = new std::vector<GRAPH_TYPE>;
	y = new std::vector<std::vector<size_t> >();

	T (*conversion)(T) = NULL;//&convFunc;

	for(size_t k=Rank;k<DHParams->FeatureFiles.size();k+=ClusterSize) {
		std::string* f = &DHParams->FeatureFiles[k];
		GRAPH_TYPE tmp;
		if(DHParams->ReadBinary==1) {
			ReadRegionFileBinary(f->c_str(), tmp, true, conversion);
		} else {
			ReadRegionFile(f->c_str(), tmp, true, conversion);
		}
		phi->x->push_back(tmp);

		std::vector<size_t> yTmp;
		ReadObservationFile(ReplaceExtension(*f, std::string(".observation")).c_str(), yTmp);
		y->push_back(yTmp);
	}

	for(size_t x=0,x_e=phi->x->size();x!=x_e;++x) {
		GRAPH_TYPE* ptr = &phi->x->at(x);
		for(typename GRAPH_TYPE::iterator rit=ptr->begin(),rit_e=ptr->end();rit!=rit_e;++rit) {
			typename S::MPNode* r = *rit;
			size_t maxEl = *std::max_element(r->featureIDs.begin(), r->featureIDs.end());
			NumFeatures = (maxEl>NumFeatures) ? maxEl : NumFeatures;
		}
	}
	
	size_t NumFeaturesLocal = NumFeatures;
	if(sizeof(size_t)==8) {
		MPI::COMM_WORLD.Allreduce(&NumFeaturesLocal, &NumFeatures, 1, MPI::UNSIGNED_LONG, MPI::MAX);
	} else {
		MPI::COMM_WORLD.Allreduce(&NumFeaturesLocal, &NumFeatures, 1, MPI::UNSIGNED, MPI::MAX);
	}
	++NumFeatures;

	for(size_t x=0,x_e=phi->x->size();x!=x_e;++x) {
		GRAPH_TYPE* ptr = &phi->x->at(x);
		rBP.push_back(new S(ptr, DHParams->epsilon));
		rBP.back()->Initialize();
	}
}

template <class T, class S>
libDHSampPMultM<T, S>::libDHSampPMultM(struct libDHBaseMultM<T>::SPParameters* params, std::vector<GRAPH_TYPE>* graph, std::vector<std::vector<size_t> >* obs, bool marginalize) {
	DHParams = params;

	MPI::Init(DHParams->argc, DHParams->argv);
	ClusterSize = MPI::COMM_WORLD.Get_size();
	Rank = MPI::COMM_WORLD.Get_rank();

	NumFeatures = 0;
	for(size_t x=0,x_e=graph->size();x!=x_e;++x) {
		GRAPH_TYPE* ptr = &graph->at(x);
		for(typename GRAPH_TYPE::iterator rit=ptr->begin(),rit_e=ptr->end();rit!=rit_e;++rit) {
			typename S::MPNode* r = *rit;
			size_t maxEl = *std::max_element(r->featureIDs.begin(), r->featureIDs.end());
			NumFeatures = (maxEl>NumFeatures) ? maxEl : NumFeatures;
		}
	}
	++NumFeatures;

	phi = new struct Sample<GRAPH_TYPE>;
	y = obs;
	isMyObs = false;
	if(marginalize) {
		isMyGraph = true;
		phi->x = new std::vector<GRAPH_TYPE>;
		Marginalize(graph);
	} else  {
		isMyGraph = false;
		phi->x = graph;
	}

	for(size_t x=0,x_e=phi->x->size();x!=x_e;++x) {
		GRAPH_TYPE* ptr = &phi->x->at(x);
		rBP.push_back(new S(ptr, DHParams->epsilon));
		rBP.back()->Initialize();
	}
}

template <class T, class S>
void libDHSampPMultM<T, S>::GetAncestors(std::set<typename S::MPNode*>& Ancestors, typename S::MPNode* r) {
	for(typename std::vector<typename S::MsgContainer>::iterator c=r->Parents.begin(),c_e=r->Parents.end();c!=c_e;++c) {
		typename S::MPNode* ptr = c->node;
		Ancestors.insert(ptr);
		GetAncestors(Ancestors, ptr);
	}
}

template <class T, class S>
int libDHSampPMultM<T, S>::Marginalize(std::vector<GRAPH_TYPE>* graph) {
	ObservedEmpiricalMean.assign(NumFeatures, T(0.0));

	phi->x->assign(graph->size(), GRAPH_TYPE());
	for(size_t x=0,x_e=graph->size();x!=x_e;++x) {
		std::vector<size_t>* y_x = &y->at(x);

		std::map<typename S::MPNode*,typename S::MPNode*> OriginalRegion2NewRegion;
		std::vector<std::pair<std::set<size_t>,typename S::MPNode*> > VariableIXs2Node;

		GRAPH_TYPE* ptr = &graph->at(x);
		GRAPH_TYPE* ptrH = &phi->x->at(x);
		for(typename GRAPH_TYPE::iterator rit=ptr->begin(),rit_e=ptr->end();rit!=rit_e;++rit) {
			typename S::MPNode* r = *rit;

			size_t pos = 0;
			size_t numLatent = 0;
			for(size_t v=0,v_e=r->num_variables;v!=v_e;++v) {
				size_t x_v = (*y_x)[r->var_ix[v]];
				if(x_v==size_t(-1)) {
					++numLatent;
				} else {
					pos += x_v*r->cum_num_states[v];
				}
			}

			if(numLatent==0) {
				for(size_t k=0,k_e=r->featureIDs.size();k!=k_e;++k) {
					ObservedEmpiricalMean[r->featureIDs[k]] += (*r->potOrig[k])[pos];
				}
			} else {
				std::set<size_t> varIXset;
				for(size_t v=0,v_e=r->num_variables;v!=v_e;++v) {
					size_t x_v = (*y_x)[r->var_ix[v]];
					if(x_v==size_t(-1)) {
						varIXset.insert(r->var_ix[v]);
					}
				}
				typename S::MPNode* ExistingNode = NULL;
				for(typename std::vector<std::pair<std::set<size_t>,typename S::MPNode*> >::iterator xs=VariableIXs2Node.begin(),xs_e=VariableIXs2Node.end();xs!=xs_e;++xs) {
					if(xs->first.size()==varIXset.size()) {
						std::set<size_t> res;
						std::set_difference(xs->first.begin(),xs->first.end(),varIXset.begin(),varIXset.end(),std::inserter(res,res.end()));
						if(res.size()==0) {
							ExistingNode = xs->second;
							break;
						}
					}
				}

				if(ExistingNode==NULL) {
					typename S::MPNode* tmp = new typename S::MPNode;
					OriginalRegion2NewRegion[r] = tmp;
					VariableIXs2Node.push_back(std::make_pair<std::set<size_t>,typename S::MPNode*>(varIXset,tmp));

					tmp->num_variables = numLatent;
					tmp->flag = r->flag;
					tmp->c_r = -std::numeric_limits<T>::max();
					tmp->SR = NULL;

					tmp->num_states = new size_t[numLatent];
					tmp->cum_num_states = new size_t[numLatent+1];tmp->cum_num_states[0] = 1;
					tmp->var_ix = new size_t[numLatent];
					size_t* multiplier = new size_t[numLatent];

					size_t cnt = 0;
					for(size_t v=0,v_e=r->num_variables;v!=v_e;++v) {
						size_t x_v = (*y_x)[r->var_ix[v]];
						if(x_v==size_t(-1)) {
							multiplier[cnt] = r->cum_num_states[v];
							tmp->num_states[cnt] = r->num_states[v];
							tmp->var_ix[cnt] = r->var_ix[v];
							tmp->cum_num_states[cnt+1] = tmp->cum_num_states[cnt]*r->num_states[v];
							++cnt;
						}
					}

					size_t potSize = tmp->cum_num_states[tmp->num_variables];
					tmp->pot = new typename S::Data;
					tmp->loss = NULL;

					size_t numFeatures = r->featureIDs.size();
					tmp->featureIDs = r->featureIDs;
					tmp->potOrig.assign(numFeatures, NULL);
					for(size_t r_=0;r_<numFeatures;++r_) {
						tmp->potOrig[r_] = new typename S::Data;
						T* buf = new T[potSize];
						typename S::Data* oldbuf = r->potOrig[r_];
						for(size_t x_v=0;x_v!=potSize;++x_v) {
							size_t oldPos = pos;
							for(size_t v=0,v_e=tmp->num_variables;v!=v_e;++v) {
								oldPos += ((x_v/tmp->cum_num_states[v])%tmp->num_states[v]) * multiplier[v];
							}

							buf[x_v] = (*oldbuf)[oldPos];
						}
						tmp->potOrig[r_]->ResetDataNew(buf,potSize,true,true);
					}

					delete [] multiplier;

					ptrH->push_back(tmp);
				} else {
					OriginalRegion2NewRegion[r] = ExistingNode;

					size_t* num_states = new size_t[numLatent];
					size_t* cum_num_states = new size_t[numLatent+1];cum_num_states[0] = 1;
					size_t* var_ix = new size_t[numLatent];
					size_t* multiplier = new size_t[numLatent];

					size_t cnt = 0;
					for(size_t v=0,v_e=r->num_variables;v!=v_e;++v) {
						size_t x_v = (*y_x)[r->var_ix[v]];
						if(x_v==size_t(-1)) {
							multiplier[cnt] = r->cum_num_states[v];
							num_states[cnt] = r->num_states[v];
							var_ix[cnt] = r->var_ix[v];
							cum_num_states[cnt+1] = cum_num_states[cnt]*r->num_states[v];
							++cnt;
						}
					}

					size_t potSize = cum_num_states[numLatent];
					size_t numFeatures = r->featureIDs.size();
					for(size_t r_=0;r_<numFeatures;++r_) {
						std::vector<size_t>::iterator rFound = std::find(ExistingNode->featureIDs.begin(), ExistingNode->featureIDs.end(), r->featureIDs[r_]);

						size_t r_tmp;
						if(rFound==ExistingNode->featureIDs.end()) {
							r_tmp = ExistingNode->featureIDs.size();
							ExistingNode->featureIDs.push_back(r->featureIDs[r_]);
							typename S::Data* buf = new typename S::Data;
							/*T* bufbuf = new T[potSize];
							std::fill(bufbuf,bufbuf+potSize, T(0.0));
							buf->ResetDataNew(bufbuf,potSize,true,true);*/
							buf->ResetDataNew((T*)NULL,(size_t*)NULL,0,true);
							ExistingNode->potOrig.push_back(buf);
						} else {
							r_tmp = rFound-ExistingNode->featureIDs.begin();
						}

						T* NewDataVec = new T[potSize];
						typename S::Data* buf = ExistingNode->potOrig[r_tmp];
						typename S::Data* oldbuf = r->potOrig[r_];

						for(size_t x_v=0;x_v!=potSize;++x_v) {
							size_t oldPos = pos;
							std::map<size_t,size_t> VarStates;
							for(size_t v=0;v!=numLatent;++v) {
								size_t buf = (x_v/cum_num_states[v])%num_states[v];
								VarStates[var_ix[v]] = buf;
								oldPos += buf * multiplier[v];
							}

							size_t newPos = 0;
							for(size_t v=0;v!=numLatent;++v) {
								newPos += VarStates[ExistingNode->var_ix[v]]*ExistingNode->cum_num_states[v];
							}

							NewDataVec[newPos] = (*buf)[newPos] + (*oldbuf)[oldPos];
						}
						buf->ResetDataNew(NewDataVec,potSize,true,true);
					}

					delete [] num_states;
					delete [] cum_num_states;
					delete [] var_ix;
					delete [] multiplier;
				}
			}
		}

		for(typename GRAPH_TYPE::iterator rit=ptr->begin(),rit_e=ptr->end();rit!=rit_e;++rit) {
			typename S::MPNode* r = *rit;
			typename std::map<typename S::MPNode*,typename S::MPNode*>::iterator res = OriginalRegion2NewRegion.find(r);
			if(res!=OriginalRegion2NewRegion.end()) {
				size_t potSize = res->second->cum_num_states[res->second->num_variables];
				for(typename std::vector<typename S::MsgContainer>::iterator p=r->Parents.begin(),p_e=r->Parents.end();p!=p_e;++p) {
					typename std::map<typename S::MPNode*,typename S::MPNode*>::iterator res_p = OriginalRegion2NewRegion.find(p->node);
					assert(res_p!=OriginalRegion2NewRegion.end());

					if(res_p->second!=res->second) {
						typename S::MsgContainer ptmp;
						ptmp.node = res_p->second;
						ptmp.lambda = new T[potSize];
						memset(ptmp.lambda, 0, sizeof(T)*potSize);
						res->second->Parents.push_back(ptmp);

						typename S::MsgContainer ctmp;
						ctmp.node = res->second;
						ctmp.lambda = ptmp.lambda;
						res_p->second->Childs.push_back(ctmp);
					}
				}
			}
		}

		if(DHParams->BetheCountingNumbers>0) {
			bool AllRegionsAssigned = false;
			while(!AllRegionsAssigned) {
				AllRegionsAssigned = true;
				for(typename GRAPH_TYPE::iterator rit=ptrH->begin(),rit_e=ptrH->end();rit!=rit_e;++rit) {
					typename S::MPNode* r = *rit;

					if(r->c_r!=-std::numeric_limits<T>::max()) {
						continue;
					}

					std::set<typename S::MPNode*> Ancestors;
					GetAncestors(Ancestors, r);

					bool allAncestorsAssigned = true;
					T ancestorSum = T(0.0);
					for(typename std::set<typename S::MPNode*>::iterator ac=Ancestors.begin(),ac_e=Ancestors.end();ac!=ac_e;++ac) {
						typename S::MPNode* ptr = *ac;
						if(ptr->c_r==-std::numeric_limits<T>::max()) {
							allAncestorsAssigned = false;
							break;
						}
						ancestorSum += ptr->c_r;
					}

					if(allAncestorsAssigned) {
						r->c_r = T(1.0) - ancestorSum;
					} else {
						AllRegionsAssigned = false;
					}
				}
			}
		} else {
			for(typename GRAPH_TYPE::iterator rit=ptrH->begin(),rit_e=ptrH->end();rit!=rit_e;++rit) {
				typename S::MPNode* r = *rit;
				r->c_r = T(1.0);
			}
		}
	}
	return 0;
}

template <class T, class S>
libDHSampPMultM<T, S>::~libDHSampPMultM() {
	DHClear();
	MPI::Finalize();
}

template <class T, class S>
int libDHSampPMultM<T, S>::DHClear() {
	for(size_t x=0,x_e=rBP.size();x!=x_e;++x) {
		delete rBP[x];
	}
	rBP.clear();
	
	if(isMyObs) {
		if(y!=NULL) {
			y->clear();
			delete y;
			y = NULL;
		}
	}
	if(isMyGraph) {
		if(phi!=NULL) {
			for(size_t x=0;x<phi->x->size();++x) {
				ClearFeature(phi->x->at(x));
			}
			phi->x->clear();
			delete phi->x;
			delete phi;
			phi = NULL;
		}
	} else {
		delete phi;
		phi = NULL;
		y = NULL;
	}

	return 0;
}

template <class T, class S>
int libDHSampPMultM<T, S>::ClearFeature(std::vector<typename S::MPNode*>& Graph) {
	std::set<typename S::Data*> ToBeDeleted;
	for(typename std::vector<typename S::MPNode*>::iterator rit=Graph.begin(),rit_e=Graph.end();rit!=rit_e;++rit) {
		typename S::MPNode* r = *rit;
		if(r->cum_num_states!=NULL) {
			delete [] r->cum_num_states;
		}
		if(r->var_ix!=NULL) {
			delete [] r->var_ix;
		}
		if(r->num_states!=NULL) {
			delete [] r->num_states;
		}
		if(r->loss!=NULL) {
			ToBeDeleted.insert(r->loss);
		}
		if(r->pot!=NULL) {
			ToBeDeleted.insert(r->pot);
		}
		r->featureIDs.clear();
		for(typename std::vector<typename S::Data*>::iterator phi_r=r->potOrig.begin(),phi_r_e=r->potOrig.end();phi_r!=phi_r_e;++phi_r) {
			if(*phi_r!=NULL) {
				ToBeDeleted.insert(*phi_r);
			}
		}
		r->potOrig.clear();
		for(typename std::vector<typename S::MsgContainer>::iterator ci=r->Childs.begin(),ci_e=r->Childs.end();ci!=ci_e;++ci) {
			if(ci->lambda!=NULL) {
				delete [] ci->lambda;
			}
		}
		delete r;
	}
	for(typename std::set<typename S::Data*>::iterator iter=ToBeDeleted.begin(),iter_e=ToBeDeleted.end();iter!=iter_e;++iter) {
		delete *iter;
	}
	return 0;
}

template <class T, class S>
size_t libDHSampPMultM<T, S>::GetNumFeatures() {
	return NumFeatures;
}

template <class T, class S>
std::vector<typename libDHSampPMultM<T,S>::GRAPH_TYPE>* libDHSampPMultM<T, S>::GetGraph() {
	return phi->x;
}

template <class T, class S>
std::vector<std::vector<size_t> >* libDHSampPMultM<T, S>::GetObservation() {
	return y;
}

template <class T, class S>
int libDHSampPMultM<T, S>::ClearMessages() {
	for(size_t x=0;x<rBP.size();++x) {
		GRAPH_TYPE* ptr = rBP[x]->Graph;
		for(typename GRAPH_TYPE::iterator rit=ptr->begin(),rit_e=ptr->end();rit!=rit_e;++rit) {
			typename S::MPNode* r = *rit;
			size_t potSize = r->cum_num_states[r->num_variables];
			for(typename std::vector<typename S::MsgContainer>::iterator p=r->Parents.begin(),p_e=r->Parents.end();p!=p_e;++p) {
				memset((char*)p->lambda, 0, sizeof(T)*potSize);
			}
		}
	}
	return 0;
}

template <class T, class S>
int libDHSampPMultM<T, S>::ComputeEmpiricalMean(std::vector<T>* empMean) {
	std::vector<T> localEmpMean(empMean->size(), T(0.0));

	for(size_t x=0,x_e=phi->x->size();x!=x_e;++x) {
		std::vector<size_t>* y_x = &y->at(x);

		GRAPH_TYPE* ptr = &phi->x->at(x);
		for(typename GRAPH_TYPE::iterator rit=ptr->begin(),rit_e=ptr->end();rit!=rit_e;++rit) {
			typename S::MPNode* r = *rit;
			
			size_t pos = 0;
			for(size_t varIX=0,varIX_e=r->num_variables;varIX!=varIX_e;++varIX) {
				size_t y_var = (*y_x)[r->var_ix[varIX]];
				if(y_var!=size_t(-1)) {
					pos += r->cum_num_states[varIX]*y_var;
				}
			}

			for(size_t k=0,k_e=r->featureIDs.size();k!=k_e;++k) {
				localEmpMean[r->featureIDs[k]] += (*r->potOrig[k])[pos];
			}
		}
	}

	if(sizeof(T)==8) {
		MPI::COMM_WORLD.Allreduce(&localEmpMean[0], &((*empMean)[0]), int(NumFeatures), MPI::DOUBLE, MPI::SUM);
	} else if(sizeof(T)==4) {
		MPI::COMM_WORLD.Allreduce(&localEmpMean[0], &((*empMean)[0]), int(NumFeatures), MPI::FLOAT, MPI::SUM);
	}

	return 0;
}

template <class T, class S>
T libDHSampPMultM<T, S>::ComputeEmpiricalMean(std::vector<T>* empMean, std::vector<T>* theta) {
	empMean->assign(ObservedEmpiricalMean.begin(), ObservedEmpiricalMean.end());
	T retval =  ComputeFandGSpecial(theta, empMean, 3, DHParams->MPIterations);

	if(sizeof(T)==8) {
		MPI::COMM_WORLD.Bcast(&(*empMean)[0], int(empMean->size()), MPI::DOUBLE, 0);
	} else if(sizeof(T)==4) {
		MPI::COMM_WORLD.Bcast(&(*empMean)[0], int(empMean->size()), MPI::FLOAT, 0);
	}

	return retval;
}

template <class T, class S>
int libDHSampPMultM<T, S>::ExchangeNu() {
	return 0;
}

template <class T, class S>
int libDHSampPMultM<T, S>::ExchangeNuPerSample() {
	return 0;
}

template <class T, class S>
int libDHSampPMultM<T, S>::Predict(std::vector<T> *theta, std::vector<std::map<size_t,std::vector<T> > >& Beliefs) {
	Beliefs.assign(DHParams->FeatureFiles.size(), std::map<size_t,std::vector<T> >());

	int x_e = int(phi->x->size());

#pragma omp parallel for
	for(int x=0;x<x_e;++x) {
		GRAPH_TYPE* ptr = &phi->x->at(x);
		for(typename GRAPH_TYPE::iterator rit=ptr->begin(),rit_e=ptr->end();rit!=rit_e;++rit) {
			typename S::MPNode* r = *rit;

			r->pot->Clear();
			for(size_t k=0,k_e=r->featureIDs.size();k!=k_e;++k) {
				r->pot->MultiplyAdd(r->potOrig[k], theta->at(r->featureIDs[k]));
			}
		}

		CPrecisionTimer* tmr = NULL;
		if(DHParams->Verbosity>0 && Rank==0) {
			tmr = new CPrecisionTimer;
			tmr->Start();
		}

		rBP[x]->RunMP(DHParams->MPIterations, T(0.0), 1, DHParams->Verbosity, tmr);
		rBP[x]->GetResult(Beliefs[Rank + x*ClusterSize]);

		if(tmr!=NULL) {
			delete tmr;
		}
	}
	
	MergeBeliefsToRank0(Beliefs);
	return 0;
}

template <class T, class S>
int libDHSampPMultM<T, S>::MergeBeliefsToRank0(std::vector<std::map<size_t,std::vector<T> > >& Beliefs) {
	std::vector<std::vector<T> > vec_dat(Beliefs.size(), std::vector<T>());
	for(typename std::vector<std::map<size_t, std::vector<T> > >::iterator x=Beliefs.begin(),x_e=Beliefs.end();x!=x_e;++x) {
		if(x->size()>0) {
			std::vector<T>* el = &vec_dat[x-Beliefs.begin()];
			el->push_back( T(x->size()) );

			for(typename std::map<size_t, std::vector<T> >::iterator v=x->begin();v!=x->end();++v) {
				el->push_back( T(v->first) );
				el->push_back( T(v->second.size()) );

				for(typename std::vector<T>::iterator v_y=v->second.begin();v_y!=v->second.end();++v_y) {
					el->push_back(*v_y);
				}
			}
		}
	}

	int* SendSize = NULL;
	if(Rank==0) {
		SendSize = new int[ClusterSize];
	}
	for(size_t k=0;k<(Beliefs.size()+ClusterSize-1)/ClusterSize;++k) {
		int mySize = 0;
		if(k*ClusterSize+Rank<vec_dat.size()) {
			mySize = int(vec_dat[k*ClusterSize + Rank].size());
		}
		MPI::COMM_WORLD.Gather(&mySize, 1, MPI::INT, SendSize, 1, MPI::INT, 0);
		
		int* displs = NULL;
		T* rdbuf = NULL;
		if(Rank==0) {
			displs = new int[ClusterSize];
			displs[0] = 0;
			for(int m=1;m<ClusterSize;++m) {
				displs[m] = displs[m-1]+SendSize[m-1];
			}
			rdbuf = new T[displs[ClusterSize-1]+SendSize[ClusterSize-1]];
		}

		if(k*ClusterSize+Rank<vec_dat.size()) {
			if(sizeof(T)==8) {
				MPI::COMM_WORLD.Gatherv(&vec_dat[k*ClusterSize + Rank][0], mySize, MPI::DOUBLE, rdbuf, SendSize, displs, MPI::DOUBLE, 0);
			} else {
				MPI::COMM_WORLD.Gatherv(&vec_dat[k*ClusterSize + Rank][0], mySize, MPI::FLOAT, rdbuf, SendSize, displs, MPI::FLOAT, 0);
			}
		} else {
			if(sizeof(T)==8) {
				MPI::COMM_WORLD.Gatherv(NULL, 0, MPI::DOUBLE, rdbuf, SendSize, displs, MPI::DOUBLE, 0);
			} else {
				MPI::COMM_WORLD.Gatherv(NULL, 0, MPI::FLOAT, rdbuf, SendSize, displs, MPI::FLOAT, 0);
			}
		}

		if(Rank==0) {
			for(int m=0;m<ClusterSize;++m) {
				if(SendSize[m]==0) {
					continue;
				}
				T* tmp = rdbuf + displs[m];

				size_t num_v = size_t(*tmp++);
				
				for(size_t n=0;n<num_v;++n) {
					size_t v_id = size_t(*tmp++);
					size_t v_len = size_t(*tmp++);

					Beliefs[k*ClusterSize+m][v_id].assign(tmp, tmp+v_len);

					tmp += v_len;
				}
			}
			
			delete [] rdbuf;
			delete [] displs;
		}
	}
	if(Rank==0) {
		delete [] SendSize;
	}
	return 0;
}

template <class T, class S>
T libDHSampPMultM<T, S>::ComputeFandGSpecial(std::vector<T> *theta, std::vector<T> *grad, int ToReturn, size_t MPIterations) {
	T dual = T(0.0);
	std::vector<T>* machine_grad = NULL;
	if(grad!=NULL) {
		machine_grad = new std::vector<T>(grad->size(), T(0.0));
	}

#pragma omp parallel
{
	T thread_dual = T(0.0);
	std::vector<T>* thread_grad = NULL;
	if(grad!=NULL) {
		thread_grad = new std::vector<T>(grad->size(), T(0.0));
	}
	int x_e = int(phi->x->size());
#pragma omp for nowait
	for(int x=0;x<x_e;++x) {
		GRAPH_TYPE* ptr = &phi->x->at(x);
		for(typename GRAPH_TYPE::iterator rit=ptr->begin(),rit_e=ptr->end();rit!=rit_e;++rit) {
			typename S::MPNode* r = *rit;

			r->pot->Clear();
			if(r->loss!=NULL) {
				r->pot->MultiplyAdd(r->loss, T(1.0));
			}
			for(size_t k=0,k_e=r->featureIDs.size();k!=k_e;++k) {
				r->pot->MultiplyAdd(r->potOrig[k], theta->at(r->featureIDs[k]));
			}
		}

		if(DHParams->ReuseMessagesForF && grad==NULL && ToReturn==1) {
			thread_dual += rBP[x]->ComputeDual();
		} else {
			rBP[x]->RunMP(MPIterations, T(0.0), ToReturn, 0);
			thread_dual += rBP[x]->ComputeDual();
		}

		if(grad!=NULL) {
			std::map<size_t, std::vector<T> > Beliefs;
			rBP[x]->GetResult(Beliefs);
			for(typename GRAPH_TYPE::iterator rit=ptr->begin(),rit_e=ptr->end();rit!=rit_e;++rit) {
				typename S::MPNode* r = *rit;

				size_t potSize = r->cum_num_states[r->num_variables];
				std::vector<T>* bel = &Beliefs[r->flag];
				for(size_t k=0,k_e=r->featureIDs.size();k!=k_e;++k) {
					T sum = T(0.0);
					typename S::Data* tmp = r->potOrig[k];
					for(size_t r_x=0;r_x!=potSize;++r_x) {
						sum += (*bel)[r_x]*(*tmp)[r_x];
					}
					thread_grad->at(r->featureIDs[k]) += sum;
				}
			}
		}
	}
#pragma omp critical
{
	dual += thread_dual;
	if(grad!=NULL) {
		for(size_t r=0,r_e=grad->size();r!=r_e;++r) {
			machine_grad->at(r) += thread_grad->at(r);
		}
		delete thread_grad;
	}
}
}

	T finalPrimalValue = T(0.0);
	if(sizeof(T)==8) {
		MPI::COMM_WORLD.Reduce(&dual, &finalPrimalValue, 1, MPI::DOUBLE, MPI::SUM, 0);
	} else if(sizeof(T)==4) {
		MPI::COMM_WORLD.Reduce(&dual, &finalPrimalValue, 1, MPI::FLOAT, MPI::SUM, 0);
	}

	if(grad!=NULL) {
		if(sizeof(T)==8) {
			MPI::COMM_WORLD.Reduce(&(*machine_grad)[0], &(*grad)[0], int(grad->size()), MPI::DOUBLE, MPI::SUM, 0);
		} else if(sizeof(T)==4) {
			MPI::COMM_WORLD.Reduce(&(*machine_grad)[0], &(*grad)[0], int(grad->size()), MPI::FLOAT, MPI::SUM, 0);
		}
		delete machine_grad;
	}

	return finalPrimalValue;
}

template <class T, class S>
T libDHSampPMultM<T, S>::ComputeFandG(std::vector<T> *theta, std::vector<T> *empiricalMeans, std::vector<T> *grad, int) {
	T finalPrimalValue = ComputeFandGSpecial(theta, grad, 1, DHParams->CRFMPIterations);

	if(Rank==0) {
		for(size_t r=0;r<NumFeatures;++r) {
			finalPrimalValue += DHParams->C/DHParams->p*pow(fabs((*theta)[r]),DHParams->p);
			finalPrimalValue -= (*theta)[r]*(*empiricalMeans)[r];
			if(grad!=NULL) {
				(*grad)[r] += (DHParams->C*pow(fabs((*theta)[r]), DHParams->p-1)*sign2((*theta)[r]) - (*empiricalMeans)[r]);
			}
		}
	}

	return finalPrimalValue;
}

template <class T, class S>
int libDHSampPMultM<T, S>::ReadRegionFile(const char* fn, std::vector<typename S::MPNode*>& Graph, bool WithFeatureID, T (*conversion)(T)) {
	std::string NetType;
	
	std::ifstream ifs(fn, std::ios_base::in);
	if(!ifs.is_open()) {
		std::cout << "File not open. " << fn << std::endl;
		return -1;
	}

	ifs >> NetType;
	if(NetType.compare("MARKOV")!=0) {
		std::cout << "The Network type is not MARKOV." << std::endl;
		return -1;
	}

	size_t num_total_vars;
	ifs >> num_total_vars;

	std::vector<size_t> VarSizes(num_total_vars);
	for(size_t k=0;k<num_total_vars;++k) {
		ifs >> VarSizes[k];
	}
	std::vector<int> MachineIDs(num_total_vars);
	for(size_t k=0;k<num_total_vars;++k) {
		ifs >> MachineIDs[k];
	}

	std::vector<std::vector<size_t> > parentsIDs;

	std::map<size_t,size_t> FileIDToGraphID;
	char remainingLine[256];

	size_t id_ = size_t(-1);
	ifs >> id_;
	while(id_!=size_t(-1) && !ifs.eof()) {
		typename S::MPNode* tmp = new typename S::MPNode;
		FileIDToGraphID[id_] = Graph.size();

		ifs >> tmp->c_r;
		ifs >> tmp->num_variables;

		tmp->flag = int(id_);
		tmp->SR = NULL;

		tmp->num_states = new size_t[tmp->num_variables];
		tmp->cum_num_states = new size_t[tmp->num_variables+1];tmp->cum_num_states[0] = 1;
		tmp->var_ix = new size_t[tmp->num_variables];

		for(size_t v=0,v_e=tmp->num_variables;v!=v_e;++v) {
			size_t v_tmp;
			ifs >> v_tmp;
			tmp->var_ix[v] = v_tmp;
			tmp->num_states[v] = VarSizes[v_tmp];
			tmp->cum_num_states[v+1] = tmp->cum_num_states[v]*VarSizes[v_tmp];
		}

		size_t potSize = tmp->cum_num_states[tmp->num_variables];
		tmp->pot = new typename S::Data;

		if(WithFeatureID) {
			size_t numFeatures;
			ifs >> numFeatures;
			tmp->featureIDs.assign(numFeatures, 0);
			tmp->potOrig.assign(numFeatures, NULL);
			for(size_t r=0;r<numFeatures;++r) {
				ifs >> tmp->featureIDs[r];
				tmp->potOrig[r] = new typename S::Data;
				T* buf = new T[potSize];
				for(size_t x_v=0;x_v!=potSize;++x_v) {
					ifs >> buf[x_v];
					if(conversion!=NULL) {
						buf[x_v] = (*conversion)(buf[x_v]);
					}
				}
				tmp->potOrig[r]->ResetDataNew(buf,potSize,true,true);
			}
			tmp->loss = new typename S::Data;
			T* buf = new T[potSize];
			for(size_t x_v=0;x_v!=potSize;++x_v) {
				ifs >> buf[x_v];
				if(conversion!=NULL) {
					buf[x_v] = (*conversion)(buf[x_v]);
				}
			}
			tmp->loss->ResetDataNew(buf,potSize,true,true);
		} else {
			tmp->loss = NULL;
			T* buf = new T[potSize];
			for(size_t x_v=0;x_v!=potSize;++x_v) {
				ifs >> buf[x_v];
				if(conversion!=NULL) {
					buf[x_v] = (*conversion)(buf[x_v]);
				}
			}
			tmp->pot->ResetDataNew(buf,potSize,true,true);
		}

		size_t numParents;
		ifs >> numParents;
		std::vector<size_t> ParentConnections(numParents, 0);
		for(size_t p=0;p!=numParents;++p) {
			ifs >> ParentConnections[p];
		}
		parentsIDs.push_back(ParentConnections);

		Graph.push_back(tmp);

		ifs.getline(remainingLine, 255, '\n');

		id_ = size_t(-1);
		ifs >> id_;
	}

	ifs.close();

	for(size_t k=0,k_e=parentsIDs.size();k!=k_e;++k) {
		typename S::MPNode* child_node = Graph[k];
		size_t potSize = child_node->cum_num_states[child_node->num_variables];

		for(size_t m=0,m_e=parentsIDs[k].size();m!=m_e;++m) {
			typename S::MPNode* parent_node = Graph[FileIDToGraphID[parentsIDs[k][m]]];

			typename S::MsgContainer ptmp;
			ptmp.node = parent_node;
			ptmp.lambda = new T[potSize];
			memset(ptmp.lambda, 0, sizeof(T)*potSize);
			child_node->Parents.push_back(ptmp);

			typename S::MsgContainer ctmp;
			ctmp.node = child_node;
			ctmp.lambda = ptmp.lambda;
			parent_node->Childs.push_back(ctmp);
		}
	}
	
	return 0;
}

template <class T, class S>
int libDHSampPMultM<T, S>::ReadRegionFileBinary(const char* fn, std::vector<typename S::MPNode*>& Graph, bool WithFeatureID, T (*conversion)(T)) {
	std::ifstream ifs(fn, std::ios_base::in |  std::ios_base::binary);
	if(!ifs.is_open()) {
		std::cout << "File not open. " << fn << std::endl;
		return -1;
	}

	int num_total_vars;
	ifs.read((char*)&num_total_vars, sizeof(int));

	std::vector<int> VarSizes(num_total_vars);
	ifs.read((char*)&VarSizes[0], num_total_vars*sizeof(int));
	std::vector<int> MachineIDs(num_total_vars);
	ifs.read((char*)&MachineIDs[0], num_total_vars*sizeof(int));

	std::vector<std::vector<int> > parentsIDs;

	std::map<int,size_t> FileIDToGraphID;
	std::map<int,std::pair<typename S::Data*,size_t> > FeatureIDToDataPtr;

	int num_regions;
	ifs.read((char*)&num_regions, sizeof(int));
	for(int k=0;k<num_regions;++k) {
		typename S::MPNode* tmp = new typename S::MPNode;

		int id_;
		ifs.read((char*)&id_, sizeof(int));
		FileIDToGraphID[id_] = Graph.size();
		tmp->flag = id_;

		double c_r;
		ifs.read((char*)&c_r, sizeof(double));
		tmp->c_r = T(c_r);

		int num_vars;
		ifs.read((char*)&num_vars, sizeof(int));
		tmp->num_variables = size_t(num_vars);

		tmp->SR = NULL;

		tmp->num_states = new size_t[tmp->num_variables];
		tmp->cum_num_states = new size_t[tmp->num_variables+1];tmp->cum_num_states[0] = 1;
		tmp->var_ix = new size_t[tmp->num_variables];

		std::vector<int> VarIXs(num_vars);
		ifs.read((char*)&VarIXs[0], num_vars*sizeof(int));
		for(size_t v=0,v_e=tmp->num_variables;v!=v_e;++v) {
			tmp->var_ix[v] = VarIXs[v];
			tmp->num_states[v] = VarSizes[VarIXs[v]];
			tmp->cum_num_states[v+1] = tmp->cum_num_states[v]*tmp->num_states[v];
		}

		size_t potSize = tmp->cum_num_states[tmp->num_variables];
		tmp->pot = new typename S::Data;

		if(WithFeatureID) {
			int numFeatures;
			ifs.read((char*)&numFeatures, sizeof(int));

			tmp->featureIDs.assign(numFeatures, 0);
			tmp->potOrig.assign(numFeatures, NULL);
			for(int r=0;r<numFeatures;++r) {
				int weightIX;
				ifs.read((char*)&weightIX, sizeof(int));
				tmp->featureIDs[r] = size_t(weightIX);

				int featureID;
				ifs.read((char*)&featureID, sizeof(int));

				typename std::map<int,std::pair<typename S::Data*,size_t> >::iterator DataPtr = FeatureIDToDataPtr.find(featureID);
				if(DataPtr==FeatureIDToDataPtr.end()) {
					typename S::Data* bufbuf = new typename S::Data;
					tmp->potOrig[r] = bufbuf;
					FeatureIDToDataPtr[featureID] = std::make_pair<typename S::Data*,size_t>(bufbuf,potSize);
				} else {
					tmp->potOrig[r] = DataPtr->second.first;
				}
			}

			int lossID;
			ifs.read((char*)&lossID, sizeof(int));

			typename std::map<int,std::pair<typename S::Data*,size_t> >::iterator LossPtr = FeatureIDToDataPtr.find(lossID);
			if(LossPtr==FeatureIDToDataPtr.end()) {
				typename S::Data* bufbuf = new typename S::Data;
				tmp->loss = bufbuf;
				FeatureIDToDataPtr[lossID] = std::make_pair<typename S::Data*,size_t>(bufbuf,potSize);
			} else {
				tmp->loss = LossPtr->second.first;
			}
		} else {
			tmp->loss = NULL;

			int PotentialID;
			ifs.read((char*)&PotentialID, sizeof(int));

			typename std::map<int,std::pair<typename S::Data*,size_t> >::iterator DataPtr = FeatureIDToDataPtr.find(PotentialID);
			if(DataPtr==FeatureIDToDataPtr.end()) {
				FeatureIDToDataPtr[PotentialID] = std::make_pair<typename S::Data*, size_t>(tmp->pot, potSize);
			}
		}

		int numParents;
		ifs.read((char*)&numParents, sizeof(int));
		if(numParents>0) {
			std::vector<int> ParentConnections(numParents, 0);
			ifs.read((char*)&ParentConnections[0], numParents*sizeof(int));
			parentsIDs.push_back(ParentConnections);
		}

		Graph.push_back(tmp);
	}

	int numPotentials;
	ifs.read((char*)&numPotentials, sizeof(int));
	for(int k=0;k<numPotentials;++k) {
		int potInfo[2];
		ifs.read((char*)potInfo, 2*sizeof(int));
		
		typename std::map<int,std::pair<typename S::Data*,size_t> >::iterator DataPtr = FeatureIDToDataPtr.find(potInfo[0]);
		if(DataPtr==FeatureIDToDataPtr.end()) {
			ifs.seekg(potInfo[1]*sizeof(double), std::ios_base::cur);
		} else {
			if(DataPtr->second.second!=size_t(potInfo[1])) {
				std::cout << "Potentials don't match." << std::endl;
				return -1;
			}
			double* pot = new double[size_t(potInfo[1])];
			ifs.read((char*)pot, sizeof(double)*potInfo[1]);
			T* buf = new T[potInfo[1]];
			for(int sz=0;sz<potInfo[1];++sz) {
				buf[sz] = T(pot[sz]);
				if(conversion!=NULL) {
					buf[sz] = (*conversion)(buf[sz]);
				}
			}
			DataPtr->second.first->ResetDataNew(buf,size_t(potInfo[1]),true,true);
			delete [] pot;
		}
	}

	ifs.close();

	for(size_t k=0,k_e=parentsIDs.size();k!=k_e;++k) {
		typename S::MPNode* child_node = Graph[k];
		size_t potSize = child_node->cum_num_states[child_node->num_variables];

		for(size_t m=0,m_e=parentsIDs[k].size();m!=m_e;++m) {
			typename S::MPNode* parent_node = Graph[FileIDToGraphID[parentsIDs[k][m]]];

			typename S::MsgContainer ptmp;
			ptmp.node = parent_node;
			ptmp.lambda = new T[potSize];
			memset(ptmp.lambda, 0, sizeof(T)*potSize);
			child_node->Parents.push_back(ptmp);

			typename S::MsgContainer ctmp;
			ctmp.node = child_node;
			ctmp.lambda = ptmp.lambda;
			parent_node->Childs.push_back(ctmp);
		}
	}
	
	return 0;
}

template <class T, class S>
std::string libDHSampPMultM<T, S>::ReplaceExtension(std::string orig, std::string extension) {
	std::string tmp(orig);
	tmp = orig.substr(0, orig.find_last_of(".")).append(extension);
	return tmp;
}

template <class T, class S>
int libDHSampPMultM<T, S>::ReadObservationFile(const char* fn, std::vector<size_t>& y) {
	std::ifstream ifs(fn, std::ios_base::in);
	if(!ifs.is_open()) {
		std::cout << "File not open. " << fn << std::endl;
		return -1;
	}

	size_t num_total_vars;
	ifs >> num_total_vars;

	y.assign(num_total_vars, size_t(-1));
	for(size_t k=0;k<num_total_vars;++k) {
		ifs >> y[k];
	}
	ifs.close();

	return 0;
}

template <class T, class S>
int libDHSampPMultM<T, S>::GetClusterSize() {
	return ClusterSize;
}

template <class T, class S>
int libDHSampPMultM<T, S>::GetRank() {
	return Rank;
}

template class libDHSampPMultM<double, libRegionBPP<double,false> >;
template class libDHSampPMultM<float, libRegionBPP<float,false> >;
