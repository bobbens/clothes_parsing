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
#include <map>
#include <iostream>
#include <cmath>
#include <set>
#include <stack>
#include <queue>
#include <string.h>
#include <omp.h>

#include <stdio.h>
#include <stdlib.h>

#include "libRegionBPP.h"

#define WITH_AGREEMENT

template <class T, bool B>
libRegionBPP<T,B>::libRegionBPP(std::vector<MPNode*>* Graph, T epsilon) {
	factorThreshold = T(-1e-5);
	agreementThreshold = T(1e-5);
	this->epsilon = epsilon;
	this->Graph = Graph;
	isMyGraph = false;

	//TreeDetection();

	GraphColoring(i2t<B>());
}

template <class T, bool B>
libRegionBPP<T,B>::libRegionBPP(std::map<size_t,std::pair<struct PotentialStruct,std::set<size_t> > > &LocalBPPot, std::map<size_t,std::pair<struct PotentialStruct,std::set<size_t> > > &FactorBPPot, T epsilon) {
	Graph = new std::vector<MPNode*>;
	ConvertToGraph(LocalBPPot, FactorBPPot);

	factorThreshold = T(-1e-5);
	agreementThreshold = T(1e-5);
	this->epsilon = epsilon;
	isMyGraph = true;

	//TreeDetection();
	
	GraphColoring(i2t<B>());
}

template <class T, bool B>
libRegionBPP<T,B>::~libRegionBPP() {
	if(isMyGraph) {
		for(typename std::vector<MPNode*>::iterator r=Graph->begin(),r_e=Graph->end();r!=r_e;++r) {
			if((*r)->cum_num_states!=NULL) {
				delete [] (*r)->cum_num_states;
			}
			if((*r)->var_ix!=NULL) {
				delete [] (*r)->var_ix;
			}
			if((*r)->pot!=NULL) {
				delete [] (*r)->pot;
			}
			for(typename std::vector<MsgContainer>::iterator ci=(*r)->Childs.begin(),ci_e=(*r)->Childs.end();ci!=ci_e;++ci) {
				if(ci->lambda!=NULL) {
					delete [] ci->lambda;
				}
			}
			delete *r;
		}
		if(Graph!=NULL) {
			delete Graph;
			Graph = NULL;
		}
	}

}

template <class T, bool B>
int libRegionBPP<T,B>::TreeDetection() {
	std::queue<T> Storage;
	for(typename std::vector<MPNode*>::iterator r=Graph->begin(),r_e=Graph->end();r!=r_e;++r) {
		(*r)->Color = size_t(-1);
		for(typename std::vector<MsgContainer>::iterator p_hat=(*r)->Parents.begin(),p_hat_e=(*r)->Parents.end();p_hat!=p_hat_e;++p_hat) {
			Storage.push(p_hat->lambda[0]);
			p_hat->lambda[0] = T(0.0);
		}
		for(typename std::vector<MsgContainer>::iterator c_hat=(*r)->Childs.begin(),c_hat_e=(*r)->Childs.end();c_hat!=c_hat_e;++c_hat) {
			Storage.push(c_hat->lambda[0]);
			c_hat->lambda[0] = T(0.0);
		}
	}

	std::set<MPNode*> AllNodes(Graph->begin(), Graph->end());
	std::stack<MPNode*> S;
	bool NotATree = false;
	while(AllNodes.size()>0 && NotATree==false) {
		S.push(*AllNodes.begin());
		while(S.size()>0 && NotATree==false) {
			MPNode* tmp = S.top();
			S.pop();
			for(typename std::vector<MsgContainer>::iterator p_hat=tmp->Parents.begin(),p_hat_e=tmp->Parents.end();p_hat!=p_hat_e && NotATree==false;++p_hat) {
				if(p_hat->lambda[0]==T(1.0)) {
					continue;
				}
				p_hat->lambda[0] = T(1.0);
				MPNode* ptr = p_hat->node;
				if(ptr->Color==size_t(-1)) {
					ptr->Color = 0;
					S.push(ptr);
				} else {
					NotATree = true;
				}
			}

			for(typename std::vector<MsgContainer>::iterator c_hat=tmp->Childs.begin(),c_hat_e=tmp->Childs.end();c_hat!=c_hat_e && NotATree==false;++c_hat) {
				if(c_hat->lambda[0]==T(1.0)) {
					continue;
				}
				c_hat->lambda[0] = T(1.0);
				MPNode* ptr = c_hat->node;
				if(ptr->Color==size_t(-1)) {
					ptr->Color = 0;
					S.push(ptr);
				} else {
					NotATree = true;
					break;
				}
			}
			AllNodes.erase(tmp);
		}
	}
	for(typename std::vector<MPNode*>::iterator r=Graph->begin(),r_e=Graph->end();r!=r_e;++r) {
		for(typename std::vector<MsgContainer>::iterator p_hat=(*r)->Parents.begin(),p_hat_e=(*r)->Parents.end();p_hat!=p_hat_e;++p_hat) {
			p_hat->lambda[0] = Storage.front();
			Storage.pop();
		}
		for(typename std::vector<MsgContainer>::iterator c_hat=(*r)->Childs.begin(),c_hat_e=(*r)->Childs.end();c_hat!=c_hat_e;++c_hat) {
			c_hat->lambda[0] = Storage.front();
			Storage.pop();
		}
	}
	if(NotATree) {
		ProcessingOrder.clear();
	}
	return 0;
}

template <class T, bool B>
int libRegionBPP<T,B>::GraphColoring(i2t<false>) {
	return 0;
}

template <class T, bool B>
int libRegionBPP<T,B>::GraphColoring(i2t<true>) {
	for(typename std::vector<MPNode*>::iterator r=Graph->begin(),r_e=Graph->end();r!=r_e;++r) {
		(*r)->Color = size_t(-1);
	}

	for(typename std::vector<MPNode*>::iterator rit=Graph->begin(),rit_e=Graph->end();rit!=rit_e;++rit) {
		MPNode* r = *rit;
		if(r->Parents.size()==0) {
			continue;
		}

		std::set<size_t> NeighboringColors;
		for(typename std::vector<MsgContainer>::iterator c=r->Childs.begin(),c_e=r->Childs.end();c!=c_e;++c) {
			size_t col = c->node->Color;
			if(col != size_t(-1)) {
				NeighboringColors.insert(col);
			}
		}
		for(typename std::vector<MsgContainer>::iterator p=r->Parents.begin(),p_e=r->Parents.end();p!=p_e;++p) {
			size_t col = p->node->Color;
			if(col != size_t(-1)) {
				NeighboringColors.insert(col);
			}
			MPNode* p_ptr = p->node;
			for(typename std::vector<MsgContainer>::iterator c_hat=p_ptr->Childs.begin(),c_hat_e=p_ptr->Childs.end();c_hat!=c_hat_e;++c_hat) {
				size_t col = c_hat->node->Color;
				if(col != size_t(-1)) {
					NeighboringColors.insert(col);
				}
			}
		}

		if(ColorToNodes.size()>0) {
			size_t cnt = 0;
			bool ColorAssigned = false;
			for(std::set<size_t>::iterator col_it=NeighboringColors.begin(),col_it_e=NeighboringColors.end();col_it!=col_it_e;++col_it, ++cnt) {
				if(*col_it>cnt) {
					r->Color = cnt;
					ColorToNodes[cnt].push_back(r);
					ColorAssigned = true;
					break;
				}
			}
			if(!ColorAssigned) {
				size_t sz = ColorToNodes.size();
				if(cnt<sz) {
					r->Color = cnt;
					ColorToNodes[cnt].push_back(r);
				} else {
					r->Color = sz;
					ColorToNodes.push_back(std::vector<MPNode*>(1, r));
				}
			}
		} else {
			r->Color = 0;
			ColorToNodes.push_back(std::vector<MPNode*>(1, r));
		}
	}
	return 0;
}

template <class T, bool B>
int libRegionBPP<T,B>::Initialize() {
	
	return 0;
}

template <class T, bool B>
int libRegionBPP<T,B>::ConvertToGraph(std::map<size_t, std::pair<struct PotentialStruct, std::set< size_t > > >& LocalBPPot, std::map<size_t, std::pair<struct PotentialStruct, std::set< size_t > > >& FactorBPPot) {
	std::map<size_t,size_t> OldIXToNewIX;
	for(typename std::map<size_t, std::pair<struct PotentialStruct, std::set< size_t > > >::iterator v=LocalBPPot.begin(),v_e=LocalBPPot.end();v!=v_e;++v) {
		MPNode* tmp = new struct MPNode;
		tmp->SR = NULL;
		tmp->c_r = T(1.0) - T(v->second.second.size());	//T(1.0);
		tmp->num_variables = v->second.first.num_variables;
		tmp->num_states = v->second.first.num_states;

		size_t numEL = tmp->num_variables+1;
		size_t *ptr = tmp->cum_num_states = new size_t[numEL];
		size_t *sz_Ptr = tmp->num_states;
		*ptr++ = 1;
		for(size_t *ptr_e=ptr+numEL-1;ptr!=ptr_e;++ptr) {
			*ptr = *(ptr-1)**sz_Ptr++;
		}

		tmp->pot = new Data;
		tmp->pot->ResetDataNoDataModify(v->second.first.pot, tmp->cum_num_states[numEL-1]);
		tmp->var_ix = new size_t[1];
		tmp->var_ix[0] = v->first;
		OldIXToNewIX.insert(std::make_pair<size_t,size_t>(v->first,Graph->size()));
		Graph->push_back(tmp);
	}
	for(typename std::map<size_t, std::pair<struct PotentialStruct, std::set< size_t > > >::iterator alpha=FactorBPPot.begin(),alpha_e=FactorBPPot.end();alpha!=alpha_e;++alpha) {
		MPNode* tmp = new struct MPNode;
		tmp->SR = NULL;
		tmp->c_r = T(1.0);
		tmp->num_variables = alpha->second.first.num_variables;
		tmp->num_states = alpha->second.first.num_states;

		size_t numEL = tmp->num_variables+1;
		size_t *ptr = tmp->cum_num_states = new size_t[numEL];
		size_t *sz_Ptr = tmp->num_states;
		*ptr++ = 1;
		for(size_t *ptr_e=ptr+numEL-1;ptr!=ptr_e;++ptr) {
			*ptr = *(ptr-1)**sz_Ptr++;
		}

		tmp->pot = new Data;
		tmp->pot->ResetDataNoDataModify(alpha->second.first.pot, tmp->cum_num_states[numEL-1]);
		tmp->flag = int(alpha->first);
		tmp->var_ix = new size_t[tmp->num_variables];
		ptr = tmp->var_ix;
		for(std::set<size_t>::iterator v=alpha->second.second.begin(),v_e=alpha->second.second.end();v!=v_e;++v, ++ptr) {
			*ptr = *v;

			size_t oldIX = OldIXToNewIX.find(*v)->second;
			MPNode *cnode = (*Graph)[oldIX];
			size_t numEL = cnode->num_states[0];

			MsgContainer ctmp;
			ctmp.node = cnode;
			ctmp.lambda = new T[numEL];
			memset(ctmp.lambda, 0, sizeof(T)*numEL);
			tmp->Childs.push_back(ctmp);

			MsgContainer ptmp;
			ptmp.node = tmp;
			ptmp.lambda = ctmp.lambda;
			cnode->Parents.push_back(ptmp);
		}
		Graph->push_back(tmp);
	}
	return 0;
}

template <class T, bool B>
T libRegionBPP<T,B>::Entropy(T* mem, size_t numEL) {
	T retVal = T(0.0);
	for(size_t x_r=0;x_r!=numEL;++x_r) {
		retVal += mem[x_r]*log(mem[x_r]+T(1e-20));
	}
	return -retVal;
}

template <class T, bool B>
T libRegionBPP<T,B>::ComputeMuPToR(MPNode *p, MPNode *r, std::map<size_t,size_t>& x_r) {
	size_t numVars = p->num_variables;
	size_t x_p_stat = 0;
	std::vector<size_t> newVars;
	for(size_t k=0;k<numVars;++k) {
		std::map<size_t,size_t>::iterator it=x_r.find(p->var_ix[k]);
		if(it==x_r.end()) {
			newVars.push_back(k);
		} else {
			x_p_stat += it->second*p->cum_num_states[k];
		}
	}

	size_t numEL = newVars.size();
	size_t* cumNewVars = new size_t[numEL+1];
	size_t* NewVarsSZ = new size_t[numEL];
	cumNewVars[0] = 1;
	for(std::vector<size_t>::iterator v=newVars.begin(),v_e=newVars.end(),v_b=v;v!=v_e;++v) {
		size_t pos = *v;
		size_t curPos = v-v_b;
		NewVarsSZ[curPos] = p->num_states[pos];
		cumNewVars[curPos+1] = cumNewVars[curPos]*NewVarsSZ[curPos];
	}

	T maxval = T(0.0);
	T* mem = NULL;
	T ecr = epsilon*p->c_r;
	if(ecr!=T(0.0)) {
		mem = new T[cumNewVars[numEL]];
	}
	
	for(size_t x_p=0,x_p_e=cumNewVars[numEL];x_p!=x_p_e;++x_p) {
		//individual vars;
		std::map<size_t,size_t> var(x_r);
		size_t x_p_real = x_p_stat;
		for(size_t varIX=0;varIX<newVars.size();++varIX) {
			size_t varVal = (x_p/cumNewVars[varIX])%NewVarsSZ[varIX];
			var[p->var_ix[newVars[varIX]]] = varVal;
			x_p_real += varVal*p->cum_num_states[newVars[varIX]];
		}

		T buf = (p->pot==NULL) ? T(0.0) : (*p->pot)[x_p_real];

		buf += (p->SR==NULL) ? T(0.0) : p->SR->nu[x_p_real];

		for(typename std::vector<MsgContainer>::iterator p_hat=p->Parents.begin(),p_hat_e=p->Parents.end();p_hat!=p_hat_e;++p_hat) {
			buf -= p_hat->lambda[x_p_real];
		}

		for(typename std::vector<MsgContainer>::iterator c_hat=p->Childs.begin(),c_hat_e=p->Childs.end();c_hat!=c_hat_e;++c_hat) {
			if(c_hat->node!=r) {
				size_t x_c_hat = 0;
				for(size_t k=0,k_e=c_hat->node->num_variables;k!=k_e;++k) {
					x_c_hat += var[c_hat->node->var_ix[k]]*c_hat->node->cum_num_states[k];
				}
				buf += c_hat->lambda[x_c_hat];
			}
		}

		if(ecr!=T(0.0)) {
			buf /= ecr;
			mem[x_p] = buf;
		}
		maxval = (buf>maxval || x_p==0) ? buf : maxval;
	}

	if(ecr!=T(0.0)) {
		T sumVal = exp( mem[0] - maxval );
		for(size_t x_p=1,x_p_e=cumNewVars[numEL];x_p!=x_p_e;++x_p) {
			sumVal += exp( mem[x_p] - maxval );
		}
		maxval *= ecr;
		maxval += ecr*log(sumVal);
		delete [] mem;
	}

	delete [] cumNewVars;
	delete [] NewVarsSZ;

	return maxval;
}

template <class T, bool B>
void libRegionBPP<T,B>::ComputeLambdaRtoP(MPNode* r) {
	if(r->Parents.size()==0) {
		return;
	}

	size_t ParentLocalIX;
	T* mu_p_r = new T[r->Parents.size()];
	//T* sum_lambda = new T[r->Parents.size()];
	//memset(sum_lambda, 0, sizeof(T)*r->Parents.size());

	T sum_c_p = r->c_r;
	for(typename std::vector<MsgContainer>::iterator p=r->Parents.begin(),p_e=r->Parents.end();p!=p_e;++p) {
		sum_c_p += p->node->c_r;
	}

	size_t numVars = r->num_variables;
	size_t x_r_e=r->cum_num_states[numVars];
	for(size_t x_r=0;x_r!=x_r_e;++x_r) {
		T phi_r_x_r_prime = (r->pot == NULL) ? T(0.0) : (*r->pot)[x_r];

		phi_r_x_r_prime += (r->SR == NULL) ? T(0.0) : r->SR->nu[x_r];

		//individual vars;
		std::map<size_t,size_t> var;
		for(size_t varIX=0;varIX<numVars;++varIX) {
			var[r->var_ix[varIX]] = (x_r/r->cum_num_states[varIX])%r->num_states[varIX];
		}

		for(typename std::vector<MsgContainer>::iterator c=r->Childs.begin(),c_e=r->Childs.end();c!=c_e;++c) {
			MPNode* ptr = c->node;
			size_t x_c = 0;
			for(size_t varIX=0,varIX_e=ptr->num_variables;varIX!=varIX_e;++varIX) {
				x_c += ptr->cum_num_states[varIX]*var[ptr->var_ix[varIX]];
			}
			phi_r_x_r_prime += c->lambda[x_c];
		}

		ParentLocalIX = 0;
		for(typename std::vector<MsgContainer>::iterator p=r->Parents.begin(),p_e=r->Parents.end();p!=p_e;++p, ++ParentLocalIX) {
			mu_p_r[ParentLocalIX] = ComputeMuPToR(p->node, r, var);
			phi_r_x_r_prime += mu_p_r[ParentLocalIX];
		}

		phi_r_x_r_prime /= sum_c_p;

		ParentLocalIX = 0;
		for(typename std::vector<MsgContainer>::iterator p=r->Parents.begin(),p_e=r->Parents.end();p!=p_e;++p, ++ParentLocalIX) {
			MPNode* ptr = p->node;//ptr points on parent, i.e., ptr->c_r = c_p!!!
			T value = ptr->c_r*phi_r_x_r_prime - mu_p_r[ParentLocalIX];
			p->lambda[x_r] = value;
			//sum_lambda[ParentLocalIX] += value;
		}
	}

	/*for(size_t k=0,k_e=r->Parents.size();k!=k_e;++k) {
		sum_lambda[k] /= x_r_e;
	}

	ParentLocalIX = 0;
	for(typename std::vector<MsgContainer>::iterator p=r->Parents.begin(),p_e=r->Parents.end();p!=p_e;++p, ++ParentLocalIX) {
		for(size_t x_r=0;x_r!=x_r_e;++x_r) {
			p->lambda[x_r] -= sum_lambda[ParentLocalIX];
		}
	}*/

	for(typename std::vector<MsgContainer>::iterator p=r->Parents.begin(),p_e=r->Parents.end();p!=p_e;++p) {
		for(size_t x_r=x_r_e-1;x_r!=0;--x_r) {
			p->lambda[x_r] -= p->lambda[0];
		}
		p->lambda[0] = 0;
	}

	//delete [] sum_lambda;
	delete [] mu_p_r;
}

template <class T, bool B>
void libRegionBPP<T, B>::Iterate(i2t<true>) {
	for(typename std::vector<std::vector<MPNode*> >::iterator col=ColorToNodes.begin(),col_e=ColorToNodes.end();col!=col_e;++col) {
#pragma omp parallel for
		for(int k=0;k<int(col->size());++k) {
			ComputeLambdaRtoP(col->at(k));
		}
	}
}

template <class T, bool B>
void libRegionBPP<T, B>::Iterate(i2t<false>) {
	for(typename std::vector<MPNode*>::iterator r=Graph->begin(),r_e=Graph->end();r!=r_e;++r) {
		ComputeLambdaRtoP(*r);
	}
}

template <class T, bool B>
void libRegionBPP<T, B>::IterateTree() {
	for(typename std::vector<MPNode*>::iterator r=ProcessingOrder.begin(),r_e=ProcessingOrder.end();r!=r_e;++r) {
		ComputeLambdaRtoP(*r);
	}
	for(typename std::vector<MPNode*>::reverse_iterator r=ProcessingOrder.rbegin(),r_e=ProcessingOrder.rend();r!=r_e;++r) {
		ComputeLambdaRtoP(*r);
	}
}

template <class T, bool B>
T libRegionBPP<T,B>::RunMP(size_t MaxIter, T pdGapThr, int ToReturn, int verbose, CPrecisionTimer* CTmr) {
	T dual = T(0.0);T primalAgree = T(0.0);T primal = T(0.0);T prevDual = T(0.0);T EntropySum = T(0.0);
	if(Graph->size()>0) {
		for(size_t iter=0;iter<MaxIter;++iter) {

			if(ProcessingOrder.size()>0) {
				IterateTree();
			} else {
				Iterate(i2t<B>());
			}

			if(verbose>1) {
				dual = ComputeDual();
				ComputePrimalWithAgreement(&primal, &primalAgree, &EntropySum);
				if(verbose>2) {
					std::cout << "D: " << dual << " P(all): " << primal << " P(loc): " << primalAgree << " G: " << dual-primalAgree;
				}
			}
			if(CTmr!=NULL && verbose>0) {
				double t = CTmr->Stop();
				std::cout << " T: " << t << std::endl;
			} else if(verbose>2) {
				std::cout << std::endl;
			}

			if(iter>0 && pdGapThr!=T(0.0) && verbose>1) {
				//if(fabs(prevDual - dual) < pdGapThr) {
				if(fabs((prevDual - dual)/dual) < pdGapThr)  {
					if(verbose>2) {
						std::cout << "Improvement: " << prevDual-dual << "; Gap: " << pdGapThr << std::endl;
					}
					break;
				}
			}
			prevDual = dual;

			if(ProcessingOrder.size()>0) {
				break;
			}
		}
	}
	return ( (ToReturn==0) ? dual-primalAgree : ( (ToReturn==1) ? dual : ( (ToReturn==2) ? primalAgree : EntropySum ) ) );
}

template <class T, bool B>
int libRegionBPP<T,B>::OutputMessages() {
	for(typename std::vector<MPNode*>::iterator r=Graph->begin(),r_e=Graph->end();r!=r_e;++r) {
		if((*r)->Parents.size()==0) {
			continue;
		}
		size_t numVars = 1;
		std::cout << "Variables:";
		for(size_t v=0;v<(*r)->num_variables;++v) {
			numVars *= (*r)->num_states[v];
			std::cout << " " << (*r)->var_ix[v];
		}
		std::cout << std::endl;
		for(size_t x_r=0;x_r<numVars;++x_r) {
			T val = T(0.0);
			for(typename std::vector<MsgContainer>::iterator p=(*r)->Parents.begin(),p_e=(*r)->Parents.end();p!=p_e;++p) {
				val += p->lambda[x_r];
			}
			std::cout << "  x_r = " << x_r << ": " << val << std::endl;
		}
	}
	return 0;
}

template <class T, bool B>
T libRegionBPP<T,B>::ComputePhiHat(MPNode* r, size_t x_r, size_t numVars) {
	T phi_hat = (r->pot==NULL) ? T(0.0) : (*r->pot)[x_r];

	phi_hat += (r->SR==NULL) ? T(0.0) : r->SR->nu[x_r];

	if(r->Childs.size()>0) {
		//individual vars;
		std::map<size_t,size_t> var;
		for(size_t varIX=0;varIX<numVars;++varIX) {
			var[r->var_ix[varIX]] = (x_r/r->cum_num_states[varIX])%r->num_states[varIX];
		}

		for(typename std::vector<MsgContainer>::iterator c=r->Childs.begin(),c_e=r->Childs.end();c!=c_e;++c) {
			size_t x_c = 0;
			for(size_t varIX=0,varIX_e=c->node->num_variables;varIX!=varIX_e;++varIX) {
				x_c += c->node->cum_num_states[varIX]*var[c->node->var_ix[varIX]];
			}
			phi_hat += c->lambda[x_c];
		}
	}

	for(typename std::vector<MsgContainer>::iterator p=r->Parents.begin(),p_e=r->Parents.end();p!=p_e;++p) {
		phi_hat -= p->lambda[x_r];
	}
	return phi_hat;
}

template <class T, bool B>
T libRegionBPP<T,B>::ComputeDual() {
	return ComputeDual(i2t<B>());
}

template <class T, bool B>
T libRegionBPP<T,B>::ComputeDual(i2t<true>) {
	T dual = T(0.0);
#pragma omp parallel
{
	T dual_threadLocal = T(0.0);
	int k_e = int(Graph->size());
#pragma omp for nowait
	for(int k=0;k<k_e;++k) {
		dual_threadLocal += ComputeDualDo((*Graph)[k]);
	}
#pragma omp critical
{
	dual += dual_threadLocal;
}
}
	return dual;
}

template <class T, bool B>
T libRegionBPP<T,B>::ComputeDual(i2t<false>) {
	T dual = T(0.0);
	for(typename std::vector<MPNode*>::iterator rit=Graph->begin(),rit_e=Graph->end();rit!=rit_e;++rit) {
		dual += ComputeDualDo(*rit);
	}
	return dual;
}

template <class T, bool B>
T libRegionBPP<T,B>::ComputeDualDo(MPNode* r) {
	T dual = T(0.0);
	size_t numVars = r->num_variables;

	T maxval = T(0.0);
	T* mem = NULL;
	T ecr = epsilon*r->c_r;
	if(ecr!=T(0.0)) {
		mem = new T[r->cum_num_states[numVars]];
	}

	for(size_t x_r=0,x_r_e=r->cum_num_states[numVars];x_r!=x_r_e;++x_r) {
		T phi_hat = ComputePhiHat(r,x_r, numVars);

		if(ecr!=T(0.0)) {
			phi_hat /= ecr;
			mem[x_r] = phi_hat;
		}
		maxval = (phi_hat>maxval || x_r==0) ? phi_hat : maxval;
	}

	if(ecr!=T(0.0)) {
		dual += ecr*maxval;
		T sum = exp( mem[0] - maxval );
		for(size_t x_r=1,x_r_e=r->cum_num_states[numVars];x_r!=x_r_e;++x_r) {
			sum += exp( mem[x_r] - maxval );
		}
		dual += ecr*log(sum+T(1e-20));
		delete [] mem;
	} else {
		dual += maxval;
	}
	return dual;
}

template <class T, bool B>
T libRegionBPP<T,B>::ComputePrimal() {
	T primal = T(0.0);
	for(typename std::vector<MPNode*>::iterator rit=Graph->begin(),rit_e=Graph->end();rit!=rit_e;++rit) {
		MPNode* r = *rit;
		size_t numVars = r->num_variables;
		size_t numEL = r->cum_num_states[numVars];
		
		T* mem = new T[numEL];

		ComputeRegionBelief(r, mem);

		for(size_t x_r=0;x_r!=numEL;++x_r) {
			primal += mem[x_r]*(*r->pot)[x_r];
		}

		T ecr = epsilon*r->c_r;
		if(ecr!=T(0.0)) {
			primal += ecr*Entropy(mem, numEL);
		}
		
		delete [] mem;
	}
	return primal;
}

template <class T, bool B>
int libRegionBPP<T,B>::ComputePrimalWithAgreement(T *primal, T *primalAgree, T* EntropySum) {
	*primal = *primalAgree = *EntropySum = T(0.0);
	std::vector<T> pVals;
	for(typename std::vector<MPNode*>::iterator rit=Graph->begin(),rit_e=Graph->end();rit!=rit_e;++rit) {
		MPNode* r = *rit;
		size_t numVars = r->num_variables;
		size_t numEL = r->cum_num_states[numVars];

		T* mem = new T[numEL];
		ComputeRegionBelief(r, mem);
		r->bel = mem;

		*EntropySum += Entropy(mem, numEL);

		T p = T(0.0);
		if(r->pot!=NULL) {
			for(size_t x_r=0;x_r!=numEL;++x_r) {
				p += mem[x_r]*(*r->pot)[x_r];
			}
		}

		T ecr = epsilon*r->c_r;
		if(ecr!=T(0.0)) {
			p += ecr*Entropy(mem, numEL);
		}

		*primal += p;
		pVals.push_back(p);
	}

	typename std::vector<T>::iterator pit = pVals.begin();
	for(typename std::vector<MPNode*>::iterator rit=Graph->begin(),rit_e=Graph->end();rit!=rit_e;++rit,++pit) {
		MPNode* r = *rit;

		size_t numELr = r->cum_num_states[r->num_variables];
		size_t numVarsr = r->num_variables;

#ifdef WITH_AGREEMENT
		bool ConstraintOK = true;
		T diff = T(0.0);
		for(typename std::vector<MsgContainer>::iterator c=r->Childs.begin(),c_e=r->Childs.end();c!=c_e && ConstraintOK;++c) {
			MPNode* cN = c->node;
			size_t numEL = cN->cum_num_states[cN->num_variables];
			T* margBel = new T[numEL];
			memset(margBel, 0, sizeof(T)*numEL);

			for(size_t x_r=0;x_r!=numELr;++x_r) { //swap loop with previous one for performance instead of memory efficiency
				std::map<size_t,size_t> var;
				for(size_t varIX=0;varIX<numVarsr;++varIX) {
					var[r->var_ix[varIX]] = (x_r/r->cum_num_states[varIX])%r->num_states[varIX];
				}

				size_t x_c = 0;
				for(size_t k=0,k_e=cN->num_variables;k!=k_e;++k) {
					x_c += var[cN->var_ix[k]]*cN->cum_num_states[k];
				}

				margBel[x_c] += r->bel[x_r];
			}

			for(size_t x_c=0;x_c!=numEL;++x_c) {
				diff += fabs(margBel[x_c]-cN->bel[x_c]);
				if(fabs(margBel[x_c]-cN->bel[x_c])>agreementThreshold) {
					ConstraintOK = false;
					break;
				}
			}

			delete [] margBel;
		}

#else
		bool ConstraintOK = (r->Childs.size()==0);
#endif

		if(ConstraintOK) {
			*primalAgree += *pit;
		} else {
			for(size_t x_r=0;x_r!=numELr;++x_r) {
				std::map<size_t,size_t> var;
				for(size_t varIX=0;varIX<numVarsr;++varIX) {
					var[r->var_ix[varIX]] = (x_r/r->cum_num_states[varIX])%r->num_states[varIX];
				}

				r->bel[x_r] = T(1.0);//requires Graph to be hierarchically sorted, i.e., all childs of r have a lower index/position
				for(typename std::vector<MsgContainer>::iterator c=r->Childs.begin(),c_e=r->Childs.end();c!=c_e;++c) {
					MPNode* cN = c->node;
				
					size_t x_c = 0;
					for(size_t k=0,k_e=cN->num_variables;k!=k_e;++k) {
						x_c += var[cN->var_ix[k]]*cN->cum_num_states[k];
					}

					r->bel[x_r] *= cN->bel[x_c];
				}
			}

			if(r->pot!=NULL) {
				T* belPtr = r->bel;
				for(size_t k=0;k<numELr;++k) {
					*primalAgree += (*r->pot)[k]**belPtr++;
				}
				/*for(T *potPtr=r->pot->GetData(),*potPtr_e=r->pot->GetData()+numELr,*belPtr=r->bel;potPtr!=potPtr_e;) {
					*primalAgree += (*potPtr++)**belPtr++;
				}*/
			}

			T ecr = epsilon*r->c_r;
			if(ecr!=T(0.0)) {
				*primalAgree += ecr*Entropy(r->bel, numELr);
			}
		}
	}

	for(typename std::vector<MPNode*>::iterator rit=Graph->begin(),rit_e=Graph->end();rit!=rit_e;++rit) {
		delete [] (*rit)->bel;
	}
	return 0;
}

template <class T, bool B>
int libRegionBPP<T,B>::ComputeRegionBelief(MPNode *r, T *mem) {
	size_t numVars = r->num_variables;
	size_t numEL = r->cum_num_states[numVars];

	T maxval = T(0.0);
	T sum = T(0.0);
	T ecrORp = T(0.0);

	if(r->Parents.size()==0) {
		ecrORp = epsilon*r->c_r;

		for(size_t x_r=0;x_r!=numEL;++x_r) {
			mem[x_r] = ComputePhiHat(r, x_r, numVars);
			if(ecrORp!=T(0.0)) {
				mem[x_r] /= ecrORp;
			}

			maxval = (mem[x_r]>maxval || x_r==0) ? mem[x_r] : maxval;
		}
	} else {/*more robust in the non-convex setting*/
		ecrORp = epsilon*r->Parents[0].node->c_r;
		for(size_t x_r=0;x_r!=numEL;++x_r) {
			//individual vars;
			std::map<size_t,size_t> var;
			for(size_t varIX=0;varIX<numVars;++varIX) {
				var[r->var_ix[varIX]] = (x_r/r->cum_num_states[varIX])%r->num_states[varIX];
			}

			mem[x_r] = r->Parents[0].lambda[x_r] + ComputeMuPToR(r->Parents[0].node, r, var);
			if(ecrORp!=T(0.0)) {
				mem[x_r] /= ecrORp;
			}

			maxval = (mem[x_r]>maxval || x_r==0) ? mem[x_r] : maxval;
		}
	}

	for(size_t x_r=0;x_r!=numEL;++x_r) {
		if(ecrORp==T(0.0)) {
			mem[x_r] = (mem[x_r] - maxval)>factorThreshold ? T(1.0) : T(0.0);
		} else {
			mem[x_r] = exp( mem[x_r] - maxval );
		}
		sum += mem[x_r];
	}
	for(size_t x_r=0;x_r!=numEL;++x_r) {
		mem[x_r] /= sum;
	}

	return 0;
}

template <class T, bool B>
int libRegionBPP<T,B>::GetResult(std::map<size_t,std::vector<T> >& Beliefs) {
	return GetResult(i2t<B>(), Beliefs);
}

template <class T, bool B>
int libRegionBPP<T,B>::GetResult(i2t<true>, std::map<size_t,std::vector<T> >& Beliefs) {
	for(typename std::vector<MPNode*>::iterator rit=Graph->begin(),rit_e=Graph->end();rit!=rit_e;++rit) {
		MPNode* r = *rit;
		size_t numVars = r->num_variables;
		size_t numEL = r->cum_num_states[numVars];
		std::vector<T> tmp(numEL, T(0.0));
		Beliefs[r->flag] = tmp;
	}
#pragma omp parallel
{
	int k_e=int(Graph->size());
#pragma omp for nowait
	for(int k=0;k<k_e;++k) {
		MPNode* r = (*Graph)[k];
		ComputeRegionBelief(r, &Beliefs[r->flag][0]);
	}
}
	return 0;
}

template <class T, bool B>
int libRegionBPP<T,B>::GetResult(i2t<false>, std::map<size_t,std::vector<T> >& Beliefs) {
	for(typename std::vector<MPNode*>::iterator rit=Graph->begin(),rit_e=Graph->end();rit!=rit_e;++rit) {
		MPNode* r = *rit;
		size_t numVars = r->num_variables;
		size_t numEL = r->cum_num_states[numVars];

		std::vector<T> tmp(numEL, T(0.0));
		ComputeRegionBelief(r, &tmp[0]);
		Beliefs[r->flag] = tmp;
	}
	return 0;
}

template class libRegionBPP<double,true>;
template class libRegionBPP<float,true>;
template class libRegionBPP<double,false>;
template class libRegionBPP<float,false>;