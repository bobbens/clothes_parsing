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

#include <string.h>
#include <omp.h>
#include <mpi.h>

#include "libDHModlPMultM.h"

#define sign2(x) (( (x) > 0 ) - ( (x) < 0 ))
#define NOT_USED(x) if((x)){}

template <class T, class S>
libDHModlPMultM<T, S>::libDHModlPMultM(struct libDHBaseMultM<T>::SPParameters* params) {
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
	SharedRegions = new struct Sample<std::vector<typename S::SharedRegion*>*>;
	SharedRegions->x = new std::vector<std::vector<typename S::SharedRegion*>*>;

	T (*conversion)(T) = NULL;//&convFunc;
	int error = 0;

	for(size_t k=0;k<DHParams->FeatureFiles.size();++k) {
		std::string* f = &DHParams->FeatureFiles[k];
		GRAPH_TYPE tmp;
		std::vector<typename S::SharedRegion*>* SR = new std::vector<typename S::SharedRegion*>;
		int retVal = 0;
		if(DHParams->ReadBinary==1) {
			retVal = ReadRegionFileBinary(f->c_str(), tmp, *SR, true, conversion);
		} else {
			retVal = ReadRegionFile(f->c_str(), tmp, *SR, true, conversion);
		}
		if(retVal!=0) {
			error = 1;
			break;
		}
		phi->x->push_back(tmp);
		SharedRegions->x->push_back(SR);

		std::vector<size_t> yTmp;
		ReadObservationFile(ReplaceExtension(*f, std::string(".observation")).c_str(), yTmp);
		y->push_back(yTmp);
	}

	int errorAll = 0;
	MPI::COMM_WORLD.Allreduce(&error, &errorAll, 1, MPI::INT, MPI::MAX);
	if(errorAll>0) {
		DHClear();
		MPI::Finalize();
		throw errorAll;
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
libDHModlPMultM<T, S>::libDHModlPMultM(struct libDHBaseMultM<T>::SPParameters* params, std::vector<GRAPH_TYPE>* graph, std::vector<std::vector<size_t> >* obs, bool marginalize) {
	NOT_USED(params)
	NOT_USED(graph);
	NOT_USED(obs);
	NOT_USED(marginalize);
	assert(false);//not supported at this point in time
}

template <class T, class S>
libDHModlPMultM<T, S>::~libDHModlPMultM() {
	DHClear();
	MPI::Finalize();
}

template <class T, class S>
int libDHModlPMultM<T, S>::DHClear() {
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

	if(SharedRegions!=NULL) {
		for(typename std::vector<std::vector<typename S::SharedRegion*>*>::iterator x=SharedRegions->x->begin(),x_e=SharedRegions->x->end();x!=x_e;++x) {
			for(typename std::vector<typename S::SharedRegion*>::iterator iter=(*x)->begin(),iter_e=(*x)->end();iter!=iter_e;++iter) {
				delete [] (*iter)->nu;
				delete *iter;
			}
			(*x)->clear();
			delete *x;
		}
		SharedRegions->x->clear();
		delete SharedRegions->x;
		delete SharedRegions;
		SharedRegions = NULL;
	}

	return 0;
}

template <class T, class S>
int libDHModlPMultM<T, S>::ClearFeature(std::vector<typename S::MPNode*>& Graph) {
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
size_t libDHModlPMultM<T, S>::GetNumFeatures() {
	return NumFeatures;
}

template <class T, class S>
int libDHModlPMultM<T, S>::ClearMessages() {
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
	for(typename std::vector<std::vector<typename S::SharedRegion*>*>::iterator x=SharedRegions->x->begin(),x_e=SharedRegions->x->end();x!=x_e;++x) {
		for(typename std::vector<typename S::SharedRegion*>::iterator r_k=(*x)->begin(),r_k_e=(*x)->end();r_k!=r_k_e;++r_k) {
			typename S::SharedRegion* tmp = *r_k;
			memset((char*)tmp->nu, 0, sizeof(T)*tmp->cumSize);
		}
	}
	return 0;
}

template <class T, class S>
int libDHModlPMultM<T, S>::ComputeEmpiricalMean(std::vector<T>* empMean) {
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
int libDHModlPMultM<T, S>::ExchangeNuPerSample() {
	for(typename std::vector<std::vector<typename S::SharedRegion*>*>::iterator x=SharedRegions->x->begin(),x_e=SharedRegions->x->end();x!=x_e;++x) {
		size_t cnt = 0;
		for(typename std::vector<typename S::SharedRegion*>::iterator r_k=(*x)->begin(),r_k_e=(*x)->end();r_k!=r_k_e;++r_k) {
			typename S::SharedRegion* tmp = *r_k;
			cnt += tmp->cumSize;
		}

		std::vector<T> AllNu(cnt, T(0.0));
		cnt = 0;
		for(typename std::vector<typename S::SharedRegion*>::iterator r_k=(*x)->begin(),r_k_e=(*x)->end();r_k!=r_k_e;++r_k) {
			typename S::SharedRegion* tmp = *r_k;
			typename S::MPNode* r = tmp->Region;
			
			T* nu_all = &AllNu[cnt];

			if(r!=NULL) {
				size_t numVars = r->num_variables;
				for(size_t x_r=0;x_r<tmp->cumSize;++x_r) {
					for(typename std::vector<typename S::MsgContainer>::iterator p=r->Parents.begin(),p_e=r->Parents.end();p!=p_e;++p) {
						nu_all[x_r] -= p->lambda[x_r];
					}

					if(r->Childs.size()>0) {
						//individual vars;
						std::map<size_t,size_t> var;
						for(size_t varIX=0;varIX<numVars;++varIX) {
							var[r->var_ix[varIX]] = (x_r/r->cum_num_states[varIX])%r->num_states[varIX];
						}

						for(typename std::vector<typename S::MsgContainer>::iterator c=r->Childs.begin(),c_e=r->Childs.end();c!=c_e;++c) {
							size_t x_c = 0;
							for(size_t varIX=0,varIX_e=c->node->num_variables;varIX!=varIX_e;++varIX) {
								x_c += c->node->cum_num_states[varIX]*var[c->node->var_ix[varIX]];
							}
							nu_all[x_r] += c->lambda[x_c];
						}
					}

					nu_all[x_r] /= tmp->numMach;
				}
			}

			cnt += tmp->cumSize;
		}

		std::vector<T> Final(AllNu.size(), T(0.0));
		if(sizeof(T)==8) {
			MPI::COMM_WORLD.Allreduce(&AllNu[0], &Final[0], int(Final.size()), MPI::DOUBLE, MPI::SUM);
		} else if(sizeof(T)==4) {
			MPI::COMM_WORLD.Allreduce(&AllNu[0], &Final[0], int(Final.size()), MPI::FLOAT, MPI::SUM);
		}

		cnt = 0;
		for(typename std::vector<typename S::SharedRegion*>::iterator r_k=(*x)->begin(),r_k_e=(*x)->end();r_k!=r_k_e;++r_k) {
			typename S::SharedRegion* tmp = *r_k;
			typename S::MPNode* r = tmp->Region;

			if(r!=NULL) {
				memcpy(tmp->nu, &Final[cnt], tmp->cumSize*sizeof(T));

				size_t numVars = r->num_variables;
				for(size_t x_r=0;x_r<tmp->cumSize;++x_r) {
					for(typename std::vector<typename S::MsgContainer>::iterator p=r->Parents.begin(),p_e=r->Parents.end();p!=p_e;++p) {
						tmp->nu[x_r] += p->lambda[x_r];
					}

					if(r->Childs.size()>0) {
						//individual vars;
						std::map<size_t,size_t> var;
						for(size_t varIX=0;varIX<numVars;++varIX) {
							var[r->var_ix[varIX]] = (x_r/r->cum_num_states[varIX])%r->num_states[varIX];
						}

						for(typename std::vector<typename S::MsgContainer>::iterator c=r->Childs.begin(),c_e=r->Childs.end();c!=c_e;++c) {
							size_t x_c = 0;
							for(size_t varIX=0,varIX_e=c->node->num_variables;varIX!=varIX_e;++varIX) {
								x_c += c->node->cum_num_states[varIX]*var[c->node->var_ix[varIX]];
							}
							tmp->nu[x_r] -= c->lambda[x_c];
						}
					}
				}
			}

			cnt += tmp->cumSize;
		}
	}
	return 0;
}

template <class T, class S>
int libDHModlPMultM<T, S>::ExchangeNu() {
	for(typename std::vector<std::vector<typename S::SharedRegion*>*>::iterator x=SharedRegions->x->begin(),x_e=SharedRegions->x->end();x!=x_e;++x) {
		for(typename std::vector<typename S::SharedRegion*>::iterator r_k=(*x)->begin(),r_k_e=(*x)->end();r_k!=r_k_e;++r_k) {
			typename S::SharedRegion* tmp = *r_k;
			typename S::MPNode* r = tmp->Region;
			
			T* nu_all = new T[tmp->cumSize];
			memset((char*)nu_all, 0, sizeof(T)*tmp->cumSize);

			if(r!=NULL) {
				size_t numVars = r->num_variables;
				for(size_t x_r=0;x_r<tmp->cumSize;++x_r) {
					for(typename std::vector<typename S::MsgContainer>::iterator p=r->Parents.begin(),p_e=r->Parents.end();p!=p_e;++p) {
						nu_all[x_r] -= p->lambda[x_r];
					}

					if(r->Childs.size()>0) {
						//individual vars;
						std::map<size_t,size_t> var;
						for(size_t varIX=0;varIX<numVars;++varIX) {
							var[r->var_ix[varIX]] = (x_r/r->cum_num_states[varIX])%r->num_states[varIX];
						}

						for(typename std::vector<typename S::MsgContainer>::iterator c=r->Childs.begin(),c_e=r->Childs.end();c!=c_e;++c) {
							size_t x_c = 0;
							for(size_t varIX=0,varIX_e=c->node->num_variables;varIX!=varIX_e;++varIX) {
								x_c += c->node->cum_num_states[varIX]*var[c->node->var_ix[varIX]];
							}
							nu_all[x_r] += c->lambda[x_c];
						}
					}

					nu_all[x_r] /= tmp->numMach;
				}
			}

			if(sizeof(T)==8) {
				MPI::COMM_WORLD.Allreduce(nu_all, tmp->nu, int(tmp->cumSize), MPI::DOUBLE, MPI::SUM);
			} else if(sizeof(T)==4) {
				MPI::COMM_WORLD.Allreduce(nu_all, tmp->nu, int(tmp->cumSize), MPI::FLOAT, MPI::SUM);
			}

			if(r!=NULL) {
				size_t numVars = r->num_variables;
				for(size_t x_r=0;x_r<tmp->cumSize;++x_r) {
					for(typename std::vector<typename S::MsgContainer>::iterator p=r->Parents.begin(),p_e=r->Parents.end();p!=p_e;++p) {
						tmp->nu[x_r] += p->lambda[x_r];
					}

					if(r->Childs.size()>0) {
						//individual vars;
						std::map<size_t,size_t> var;
						for(size_t varIX=0;varIX<numVars;++varIX) {
							var[r->var_ix[varIX]] = (x_r/r->cum_num_states[varIX])%r->num_states[varIX];
						}

						for(typename std::vector<typename S::MsgContainer>::iterator c=r->Childs.begin(),c_e=r->Childs.end();c!=c_e;++c) {
							size_t x_c = 0;
							for(size_t varIX=0,varIX_e=c->node->num_variables;varIX!=varIX_e;++varIX) {
								x_c += c->node->cum_num_states[varIX]*var[c->node->var_ix[varIX]];
							}
							tmp->nu[x_r] -= c->lambda[x_c];
						}
					}
				}
			}
		}
	}
	return 0;
}

template <class T, class S>
int libDHModlPMultM<T, S>::Predict(std::vector<T> *theta, std::vector<std::map<size_t,std::vector<T> > >& Beliefs) {
	Beliefs.assign(DHParams->FeatureFiles.size(), std::map<size_t,std::vector<T> >());

	for(int x=0,x_e=int(phi->x->size());x<x_e;++x) {
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

		for(int outerIter=0,outerIter_e=int(DHParams->CRFOuterIterations)-1;outerIter<outerIter_e;++outerIter) {
			rBP[x]->RunMP(DHParams->MPIterations, T(0.0), 1, 0);

			if(DHParams->Verbosity>1) {
				std::vector<T> results(4, T(0.0));
				results[0] = rBP[x]->ComputeDual();
				rBP[x]->ComputePrimalWithAgreement(&results[1], &results[2], &results[3]);
				std::vector<T> AllResults(4, T(0.0));
				if(sizeof(T)==8) {
					MPI::COMM_WORLD.Reduce(&results[0], &AllResults[0], 4, MPI::DOUBLE, MPI::SUM, 0);
				} else if(sizeof(T)==4) {
					MPI::COMM_WORLD.Reduce(&results[0], &AllResults[0], 4, MPI::FLOAT, MPI::SUM, 0);
				}
				if(DHParams->Verbosity>2 && Rank==0) {
					std::cout << "D: " << AllResults[0] << " P(all): " << AllResults[1] << " P(loc): " << AllResults[2] << " G: " << AllResults[0]-AllResults[2];
				}
			}

			if(tmr!=NULL && DHParams->Verbosity>0) {
				double t = tmr->Stop();
				std::cout << " T: " << t << std::endl;
			} else if(Rank==0 && DHParams->Verbosity>2) {
				std::cout << std::endl;
			}
			ExchangeNuPerSample();
		}
		rBP[x]->RunMP(DHParams->MPIterations, T(0.0), 1, 0);
		
		if(DHParams->Verbosity>1) {
			std::vector<T> results(4, T(0.0));
			results[0] = rBP[x]->ComputeDual();
			rBP[x]->ComputePrimalWithAgreement(&results[1], &results[2], &results[3]);
			std::vector<T> AllResults(4, T(0.0));
			if(sizeof(T)==8) {
				MPI::COMM_WORLD.Reduce(&results[0], &AllResults[0], 4, MPI::DOUBLE, MPI::SUM, 0);
			} else if(sizeof(T)==4) {
				MPI::COMM_WORLD.Reduce(&results[0], &AllResults[0], 4, MPI::FLOAT, MPI::SUM, 0);
			}
			if(DHParams->Verbosity>2 && Rank==0) {
				std::cout << "D: " << AllResults[0] << " P(all): " << AllResults[1] << " P(loc): " << AllResults[2] << " G: " << AllResults[0]-AllResults[2];
			}
		}

		if(tmr!=NULL && DHParams->Verbosity>0) {
			double t = tmr->Stop();
			std::cout << " T: " << t << std::endl;
		} else if(Rank==0 && DHParams->Verbosity>2) {
			std::cout << std::endl;
		}
		
		rBP[x]->GetResult(Beliefs[x]);

		if(tmr!=NULL) {
			delete tmr;
		}
	}

	MergeBeliefsToRank0(Beliefs);
	return 0;
}

template <class T, class S>
int libDHModlPMultM<T, S>::MergeBeliefsToRank0(std::vector<std::map<size_t,std::vector<T> > >& Beliefs) {
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
	for(size_t k=0;k<Beliefs.size();++k) {
		int mySize = int(vec_dat[k].size());
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

		if(vec_dat[k].size()>0) {
			if(sizeof(T)==8) {
				MPI::COMM_WORLD.Gatherv(&vec_dat[k][0], mySize, MPI::DOUBLE, rdbuf, SendSize, displs, MPI::DOUBLE, 0);
			} else {
				MPI::COMM_WORLD.Gatherv(&vec_dat[k][0], mySize, MPI::FLOAT, rdbuf, SendSize, displs, MPI::FLOAT, 0);
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

					Beliefs[k][v_id].assign(tmp, tmp+v_len);

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
T libDHModlPMultM<T, S>::ComputeFandGSpecial(std::vector<T> *theta, std::vector<T> *grad, int ToReturn, size_t MPIterations, int doExchange) {
	T dual = T(0.0);
	std::vector<T>* machine_grad = NULL;
	if(grad!=NULL) {
		machine_grad = new std::vector<T>(grad->size(), T(0.0));
	}

	for(int x=0,x_e=int(phi->x->size());x<x_e;++x) {
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
			dual += rBP[x]->ComputeDual();
		} else {
			for(int outerIter=0,outerIter_e=int(DHParams->CRFOuterIterations)-1;outerIter<outerIter_e;++outerIter) {
				rBP[x]->RunMP(MPIterations, T(0.0), 1, 0);
				ExchangeNuPerSample();
			}
			rBP[x]->RunMP(MPIterations, T(0.0), 1, 0);
			if(doExchange>0) {
				ExchangeNuPerSample();
			}
			dual += rBP[x]->ComputeDual();
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
					machine_grad->at(r->featureIDs[k]) += sum;
				}
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
T libDHModlPMultM<T, S>::ComputeFandG(std::vector<T> *theta, std::vector<T> *empiricalMeans, std::vector<T> *grad, int flag) {
	T finalPrimalValue = ComputeFandGSpecial(theta, grad, 1, DHParams->CRFMPIterations, flag);

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
int libDHModlPMultM<T, S>::ReadRegionFile(const char* fn, std::vector<typename S::MPNode*>& Graph, std::vector<typename S::SharedRegion*>& SR, bool WithFeatureID, T (*conversion)(T)) {
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
	if(*std::max_element(MachineIDs.begin(), MachineIDs.end())>=ClusterSize) {
		std::cout << "Cluster configuration does not match problem." << std::endl;
		return -1;
	}

	std::vector<std::vector<size_t> > parentsIDs;

	std::map<size_t,size_t> FileIDToGraphID;
	char remainingLine[256];

	size_t id_ = size_t(-1);
	ifs >> id_;

	while(id_!=size_t(-1) && !ifs.eof()) {
		typename S::MPNode* tmp = new typename S::MPNode;

		ifs >> tmp->c_r;
		ifs >> tmp->num_variables;

		tmp->flag = int(id_);
		tmp->SR = NULL;

		tmp->num_states = new size_t[tmp->num_variables];
		tmp->cum_num_states = new size_t[tmp->num_variables+1];tmp->cum_num_states[0] = 1;
		tmp->var_ix = new size_t[tmp->num_variables];

		bool RegionOnThisMachine = false;
		std::set<int> ValidMachines;
		for(size_t v=0,v_e=tmp->num_variables;v!=v_e;++v) {
			size_t v_tmp;
			ifs >> v_tmp;
			tmp->var_ix[v] = v_tmp;
			tmp->num_states[v] = VarSizes[v_tmp];
			tmp->cum_num_states[v+1] = tmp->cum_num_states[v]*VarSizes[v_tmp];
			RegionOnThisMachine = RegionOnThisMachine || MachineIDs[v_tmp]==Rank;
			ValidMachines.insert(MachineIDs[v_tmp]);
		}

		if(ValidMachines.size()>1) {
			typename S::SharedRegion* sr = new typename S::SharedRegion;
			size_t potSize = tmp->cum_num_states[tmp->num_variables];
			sr->cumSize = potSize;
			sr->numMach = ValidMachines.size();
			sr->nu = new T[potSize];
			memset((char*)sr->nu, 0, potSize*sizeof(T));
			if(RegionOnThisMachine) {
				sr->Region = tmp;
				tmp->SR = sr;
			} else {
				sr->Region = NULL;
			}
			SR.push_back(sr);
		}

		if(RegionOnThisMachine) {
			FileIDToGraphID[id_] = Graph.size();

			tmp->c_r /= ValidMachines.size();

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
						buf[x_v] /= ValidMachines.size();
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
					buf[x_v] /= ValidMachines.size();
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
					buf[x_v] /= ValidMachines.size();
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
		} else {
			delete [] tmp->num_states;
			delete [] tmp->cum_num_states;
			delete [] tmp->var_ix;
			delete tmp;
		}

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
int libDHModlPMultM<T, S>::ReadRegionFileBinary(const char* fn, std::vector<typename S::MPNode*>& Graph, std::vector<typename S::SharedRegion*>& SR, bool WithFeatureID, T (*conversion)(T)) {
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

	if(*std::max_element(MachineIDs.begin(), MachineIDs.end())>=ClusterSize) {
		std::cout << "Cluster configuration does not match problem." << std::endl;
		return -1;
	}

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

		tmp->SR = NULL;

		int num_vars;
		ifs.read((char*)&num_vars, sizeof(int));
		tmp->num_variables = size_t(num_vars);

		tmp->num_states = new size_t[tmp->num_variables];
		tmp->cum_num_states = new size_t[tmp->num_variables+1];tmp->cum_num_states[0] = 1;
		tmp->var_ix = new size_t[tmp->num_variables];

		bool RegionOnThisMachine = false;
		std::set<int> ValidMachines;
		std::vector<int> VarIXs(num_vars);
		ifs.read((char*)&VarIXs[0], num_vars*sizeof(int));
		for(size_t v=0,v_e=tmp->num_variables;v!=v_e;++v) {
			tmp->var_ix[v] = VarIXs[v];
			tmp->num_states[v] = VarSizes[VarIXs[v]];
			tmp->cum_num_states[v+1] = tmp->cum_num_states[v]*tmp->num_states[v];
			RegionOnThisMachine = RegionOnThisMachine || MachineIDs[VarIXs[v]]==Rank;
			ValidMachines.insert(MachineIDs[VarIXs[v]]);
		}

		if(ValidMachines.size()>1) {
			typename S::SharedRegion* sr = new typename S::SharedRegion;
			size_t potSize = tmp->cum_num_states[tmp->num_variables];
			sr->cumSize = potSize;
			sr->numMach = ValidMachines.size();
			sr->nu = new T[potSize];
			memset((char*)sr->nu, 0, potSize*sizeof(T));
			if(RegionOnThisMachine) {
				sr->Region = tmp;
				tmp->SR = sr;
			} else {
				sr->Region = NULL;
			}
			SR.push_back(sr);
		}

		if(RegionOnThisMachine) {
			FileIDToGraphID[id_] = Graph.size();

			tmp->c_r /= ValidMachines.size();

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
		} else {
			delete [] tmp->num_states;
			delete [] tmp->cum_num_states;
			delete [] tmp->var_ix;
			delete tmp;

			if(WithFeatureID) {
				int numFeatures;
				ifs.read((char*)&numFeatures, sizeof(int));

				std::vector<int> tmp(2*numFeatures+1);
				ifs.read((char*)&tmp[0], (2*numFeatures+1)*sizeof(int));
			} else {
				int potID;
				ifs.read((char*)&potID, sizeof(int));
			}

			int numParents;
			ifs.read((char*)&numParents, sizeof(int));
			if(numParents>0) {
				std::vector<int> ParentConnections(numParents, 0);
				ifs.read((char*)&ParentConnections[0], numParents*sizeof(int));
			}
		}
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
			double* pot = new double[potInfo[1]];
			ifs.read((char*)pot, sizeof(double)*potInfo[1]);
			T* buf = new T[size_t(potInfo[1])];
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
std::string libDHModlPMultM<T, S>::ReplaceExtension(std::string orig, std::string extension) {
	std::string tmp(orig);
	tmp = orig.substr(0, orig.find_last_of(".")).append(extension);
	return tmp;
}

template <class T, class S>
int libDHModlPMultM<T, S>::ReadObservationFile(const char* fn, std::vector<size_t>& y) {
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
int libDHModlPMultM<T, S>::GetClusterSize() {
	return ClusterSize;
}

template <class T, class S>
int libDHModlPMultM<T, S>::GetRank() {
	return Rank;
}

template class libDHModlPMultM<double,libRegionBPP<double,true> >;
template class libDHModlPMultM<float,libRegionBPP<float,true> >;
