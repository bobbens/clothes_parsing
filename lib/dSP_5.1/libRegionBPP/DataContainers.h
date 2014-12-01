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
#ifndef __DATACONTAINERS_H__
#define __DATACONTAINERS_H__
#include <vector>
#include <map>
#ifdef USE_ON_WINDOWS
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif
#include <cmath>
#include <algorithm>
#include <functional>
#include <string.h>
#include <assert.h>

template <class T> class BaseClassDense;
template <class T> class BaseClassSparse;
template <class T> class BaseClassSparse2;

template <class T>
class BaseClass {
public:
	virtual ~BaseClass() {}
	virtual void Clear() = 0;
	virtual BaseClass<T>* clone() = 0;
	virtual void Set(size_t ix, T val) = 0;
	virtual T Get(size_t ix) = 0;
	virtual void Multiply(T factor) = 0;
	virtual BaseClass<T>* MultiplyAdd(BaseClass<T>*, T) = 0;
	virtual BaseClass<T>* MultiplyAdd(BaseClassDense<T>*, T) = 0;
	virtual BaseClass<T>* MultiplyAdd(BaseClassSparse<T>*, T) = 0;
	virtual BaseClass<T>* MultiplyAdd(BaseClassSparse2<T>*, T) = 0;
};

template <class T>
class BaseClassDense : public BaseClass<T> {
	friend class BaseClassSparse<T>;
	friend class BaseClassSparse2<T>;
	T* data;
	size_t sz;
	int myData;//0:not modifiable and no delete responsibility; 1: modifiable and delete responsibility
public:
	BaseClassDense(size_t numEl) {
		if(numEl>0) {
			data = new T[numEl];std::fill(data,data+numEl,T(0.0));myData = 1;sz = numEl;
		} else {
			data = NULL;myData = 1;sz = numEl;
		}
	}
	BaseClassDense(T* dataPtr, size_t numEl, int myDat) {data = dataPtr;myData = myDat;sz = numEl;}
	~BaseClassDense() {
		if(myData==1) {
			if(data!=NULL) {
				delete [] data;
				data = NULL;
			}
			sz = 0;
		}
	}
	void Clear() {
		if(myData==0 && sz>0) {
			data = new T[sz];
			memset((char*)data,0,sizeof(T)*sz);
			myData = 1;
		} else {
			if(sz>0) {
				memset((char*)data, 0, sizeof(T)*sz);
			}
		}
	}
	BaseClass<T>* clone() {
		BaseClassDense<T>* retval = new BaseClassDense<T>(sz);
		if(sz>0) {
			memcpy(retval->data, data, sizeof(T)*sz);
		}
		return retval;
	}
	void Multiply(T factor) {
		if(myData==0 && sz>0) {
			T* dataNew = new T[sz];
			memcpy((char*)dataNew,(char*)data,sizeof(T)*sz);
			myData = 1;
		}
		if(sz>0) {
			for(T *ptr=data,*ptr_e=data+sz;ptr!=ptr_e;++ptr) {
				*ptr *= factor;
			}
		}
	}
	void Set(size_t ix, T val) {
		assert(ix<sz);//sz>0
		if(myData==0) {
			T* dataNew = new T[sz];
			memcpy((char*)dataNew,(char*)data,sizeof(T)*sz);
			myData = 1;
		}
		data[ix] = val;
	}
	T Get(size_t ix) {
		assert(ix<sz);
		return data[ix];
	}
	BaseClass<T>* MultiplyAdd(BaseClass<T>* rhs, T factor) { return rhs->MultiplyAdd(this, factor);}
	BaseClass<T>* MultiplyAdd(BaseClassDense<T>* lhs, T factor) {
		assert(lhs->sz==sz && lhs->myData==1);
		for(T *rhs_data=data,*rhs_data_e=data+sz,*lhs_data=lhs->data;rhs_data!=rhs_data_e;++rhs_data,++lhs_data) {
			*lhs_data += factor**rhs_data;
		}
		return lhs;
	}
	BaseClass<T>* MultiplyAdd(BaseClassSparse<T>* lhs, T factor) {
		BaseClassDense<T>* rhs = static_cast<BaseClassDense<T>*>(clone());
		rhs->Multiply(factor);
		lhs->MultiplyAdd(rhs, T(1.0));//modifies data of rhs
		delete lhs;
		return rhs;
	}
	BaseClass<T>* MultiplyAdd(BaseClassSparse2<T>* lhs, T factor) {
		BaseClassDense<T>* rhs = static_cast<BaseClassDense<T>*>(clone());
		rhs->Multiply(factor);
		lhs->MultiplyAdd(rhs, T(1.0));//modifies data of rhs
		delete lhs;
		return rhs;
	}
};

template <class T>
class BaseClassSparse2 : public BaseClass<T> {
	friend class BaseClassDense<T>;
	friend class BaseClassSparse<T>;
	T* data;
	size_t* ix;
	size_t sz;
	int myData;//0:not modifiable;1: modifiable and delete responsibility
	void Delete() {
		if(myData==1) {
			if(data!=NULL) {delete [] data;data = NULL;}
			if(ix!=NULL) {delete [] ix;ix = NULL;}
			sz = 0;
		}
	}
public:
	BaseClassSparse2() : data(NULL), ix(NULL), sz(0), myData(1) {}
	BaseClassSparse2(T* dat, size_t* ir, size_t numEl, int myDat) : data(dat), ix(ir), sz(numEl), myData(myDat) {}
	~BaseClassSparse2() {
		Delete();
	}
	void Clear() {
		if(myData==1) {
			Delete();
		} else {
			data = NULL;
			ix = NULL;
			sz = 0;
			myData = 1;
		}
	}
	BaseClass<T>* clone() {
		if(sz!=0) {
			T* dataNew = new T[sz];
			size_t* ixNew = new size_t[sz];
			memcpy((char*)dataNew, (char*)data, sizeof(T)*sz);
			memcpy((char*)ixNew, (char*)ix, sizeof(size_t)*sz);
			return new BaseClassSparse2(dataNew, ixNew, sz, 1);
		} else {
			return new BaseClassSparse2();
		}
	}
	void Multiply(T factor) {
		if(sz>0) {
			for(T *ptr=data,*ptr_e=data+sz;ptr!=ptr_e;++ptr) {
				*ptr *= factor;
			}
		}
	}
	void Set(size_t ixEl, T val) {
		if(ix!=NULL) {
			size_t* last = ix+sz;
			size_t* el = std::lower_bound(ix,last,ixEl);//std::find_if(ix,last,std::bind2nd(std::greater_equal<T>(),ixEl));
			if(el!=last && *el==ixEl) {
				if(myData==0) {
					T* dataNew = new T[sz];
					size_t* ixNew = new size_t[sz];
					memcpy((char*)dataNew, (char*)data, sizeof(T)*sz);
					memcpy((char*)ixNew, (char*)ix, sizeof(size_t)*sz);
					data = dataNew;
					ix = ixNew;
					myData = 1;
				}
				data[el-ix] = val;
			} else {
				InsertNewBefore(el, ixEl, val);
			}
		} else {
			assert(myData==1);
			T* data = new T[1];
			size_t* ix = new size_t[1];
			data[0] = val;
			ix[0] = ixEl;
			myData = 1;
		}
	}
	T Get(size_t ixEl) {
		size_t* last = ix+sz;
		size_t* el = std::lower_bound(ix,last,ixEl);//std::find_if(ix,last,std::bind2nd(std::greater_equal<T>(),ixEl));
		if(el!=last && *el==ixEl) {
			return data[el-ix];
		} else {
			return T(0.0);
		}
	}
	BaseClass<T>* MultiplyAdd(BaseClass<T>* rhs, T factor) { return rhs->MultiplyAdd(this, factor);}
	BaseClass<T>* MultiplyAdd(BaseClassDense<T>* lhs, T factor) {
		assert(lhs->myData==1);
		T* dataptr = data;
		for(size_t *ixptr=ix,*ixptr_e=ix+sz;ixptr!=ixptr_e;++ixptr) {
			assert(*ixptr<lhs->sz);
			lhs->data[*ixptr] += factor**dataptr++;
		}
		return lhs;
	}
	BaseClass<T>* MultiplyAdd(BaseClassSparse<T>* lhs, T factor) {
		T* dataptr = data;
		for(size_t *ixptr=ix,*ixptr_e=ix+sz;ixptr!=ixptr_e;++ixptr) {
			typename BaseClassSparse<T>::CONTAINER::iterator it=lhs->data.find(*ixptr);
			if(it==lhs->data.end()) {
				lhs->data[*ixptr] = factor**dataptr++;
			} else {
				it->second += factor**dataptr++;
			}
		}
		return lhs;
	}
	BaseClass<T>* MultiplyAdd(BaseClassSparse2<T>* lhs, T factor) {
		assert(lhs->myData==1);
		if(lhs->sz>0 && sz>0) {
			T* dataptr = data;
			size_t* last = lhs->ix+lhs->sz;
			for(size_t *ixptr=ix,*ixptr_e=ix+sz;ixptr!=ixptr_e;++ixptr) {
				size_t* el = std::lower_bound(lhs->ix,last,*ixptr);//std::find_if(lhs->ix,last,std::bind2nd(std::greater_equal<T>(),*ixptr));
				if(el!=last && *el==*ixptr) {
					lhs->data[el-lhs->ix] += factor**dataptr++;
				} else {
					lhs->InsertNewBefore(el,*ixptr,factor**dataptr);
					last = lhs->ix+lhs->sz;
					++dataptr;
				}
			}
		} else if(sz>0) {
			assert(lhs->data==NULL && lhs->ix==NULL);
			lhs->data = new T[sz];
			lhs->ix = new size_t[sz];
			lhs->sz = sz;
			memcpy((char*)lhs->ix,(char*)ix,sizeof(size_t)*sz);
			for(T *dataptr=data,*dataptr_e=data+sz,*targetptr=lhs->data;dataptr!=dataptr_e;++dataptr) {
				*targetptr++ = factor**dataptr;
			}
		}
		return lhs;
	}
	void InsertNewBefore(size_t* ptr, size_t ixEntry, T valEntry) {
		size_t* ixNew = new size_t[sz+1];
		T* dataNew = new T[sz+1];
		size_t loc = ptr-ix;
		if(loc==0) {
			dataNew[0] = valEntry;
			ixNew[0] = ixEntry;
			memcpy((char*)&ixNew[1], (char*)ix, sizeof(size_t)*sz);
			memcpy((char*)&dataNew[1], (char*)data, sizeof(T)*sz);
		} else if(loc==sz) {
			memcpy((char*)ixNew, (char*)ix, sizeof(size_t)*sz);
			memcpy((char*)dataNew, (char*)data, sizeof(T)*sz);
			dataNew[sz] = valEntry;
			ixNew[sz] = ixEntry;
		} else if(loc<sz) {
			memcpy((char*)ixNew, (char*)ix, sizeof(size_t)*loc);
			memcpy((char*)dataNew, (char*)data, sizeof(T)*loc);
			ixNew[loc] = ixEntry;
			dataNew[loc] = valEntry;
			memcpy((char*)&ixNew[loc+1], (char*)&ix[loc], sizeof(size_t)*(sz-loc));
			memcpy((char*)&dataNew[loc+1], (char*)&data[loc], sizeof(T)*(sz-loc));
		} else {
			assert(false);
		}
		++sz;
		if(myData==1) {
			delete [] ix;
			delete [] data;
			ix = ixNew;
			data = dataNew;
		} else {
			ix = ixNew;
			data = dataNew;
			myData = 1;
		}
	}
};

template <class T>
class BaseClassSparse : public BaseClass<T> {
	friend class BaseClassDense<T>;
	friend class BaseClassSparse2<T>;
	typedef typename std::tr1::unordered_map<size_t,T> CONTAINER;//std::map<size_t,T>
	CONTAINER data;
public:
	BaseClassSparse() {}
	~BaseClassSparse() {}
	void Clear() {
		data.clear();
	}
	BaseClass<T>* clone() {
		BaseClassSparse<T>* retval = new BaseClassSparse<T>;
		retval->data = data;
		return retval;
	}
	void Multiply(T factor) {
		for(typename CONTAINER::iterator iter=data.begin(),iter_e=data.end();iter!=iter_e;++iter) {
			iter->second *= factor;
		}
	}
	void Set(size_t ix, T val) {
		data[ix] = val;
	}
	T Get(size_t ix) {
		typename CONTAINER::const_iterator it = data.find(ix);
		if(it==data.end()) {
			return T(0.0);
		} else {
			return it->second;
		}
	}
	BaseClass<T>* MultiplyAdd(BaseClass<T>* rhs, T factor) { return rhs->MultiplyAdd(this, factor);}
	BaseClass<T>* MultiplyAdd(BaseClassDense<T>* lhs, T factor) {
		assert(lhs->myData==1);
		for(typename CONTAINER::const_iterator iter=data.begin(),iter_e=data.end();iter!=iter_e;++iter) {
			assert(iter->first<lhs->sz);
			lhs->data[iter->first] += factor*iter->second;
		}
		return lhs;
	}
	BaseClass<T>* MultiplyAdd(BaseClassSparse<T>* lhs, T factor) {
		for(typename CONTAINER::const_iterator iter=data.begin(),iter_e=data.end();iter!=iter_e;++iter) {
			typename CONTAINER::iterator it=lhs->data.find(iter->first);
			if(it==lhs->data.end()) {
				lhs->data[iter->first] = factor*iter->second;
			} else {
				it->second += factor*iter->second;
			}
		}
		return lhs;
	}
	BaseClass<T>* MultiplyAdd(BaseClassSparse2<T>* lhs, T factor) {
		assert(lhs->myData==1);
		if(lhs->ix!=NULL) {
			for(typename CONTAINER::const_iterator iter=data.begin(),iter_e=data.end();iter!=iter_e;++iter) {
				size_t* last = lhs->ix+lhs->sz;
				size_t* el = std::lower_bound(lhs->ix,last,iter->first);//std::find_if(lhs->ix,last,std::bind2nd(std::greater_equal<T>(),iter->first));
				if(el!=last && *el==iter->first) {
					lhs->data[el-lhs->ix] += factor*iter->second;
				} else {
					lhs->InsertNewBefore(el, iter->first, factor*iter->second);
				}
			}
		} else {
			size_t sz = data.size();
			T* dataNew = new T[sz];T* datptr = dataNew;
			size_t* ixNew = new size_t[sz];size_t* ixptr = ixNew;
			for(typename CONTAINER::const_iterator iter=data.begin(),iter_e=data.end();iter!=iter_e;++iter) {
				*datptr++ = factor*iter->second;
				*ixptr++ = iter->first;
			}
			delete lhs;
			lhs = new BaseClassSparse2<T>(dataNew, ixNew, sz, 1);
		}
		return lhs;
	}
};

template <class T>
struct ProxyObject {
	BaseClass<T> *obj;
	size_t ix;
	ProxyObject<T>( BaseClass<T> *obj, size_t ix ) : obj(obj), ix(ix) {};
	/*void operator=( const T& val ) {
	obj->Set(ix,val);
	}*/
	operator T() {
		return obj->Get(ix);
	}
};

template <class T>
class DataContainer {
private:
	BaseClass<T>* obj;
public:
	DataContainer() : obj(NULL) {};
	~DataContainer() {
		Delete();
	}
	void Delete() {
		if(obj!=NULL) delete obj;
	}
	void Clear() {
		if(obj!=NULL) obj->Clear();
	}
	void MultiplyAdd(DataContainer<T>* rhs, T factor) {
		if(obj==NULL) {
			obj = rhs->obj->clone();
			obj->Multiply(factor);
		} else {
			obj = obj->MultiplyAdd(rhs->obj, factor);
		}
	}
	void ResetDataNoDataModify(T* d, size_t numEl) {//e.g., for matlab
		if(obj!=NULL) {
			delete obj;
		}
		obj = new BaseClassDense<T>(d, numEl, 0);
	}
	void ResetDataNew(T* d, size_t* ix, size_t numEl, bool takeOverDataResponsibility) {
		if(obj!=NULL) {
			delete obj;
		}
		obj = new BaseClassSparse2<T>(d,ix,numEl,((takeOverDataResponsibility)?1:0));
	}
	void ResetDataNew(T* d, size_t numEl, bool CheckForSparsity, bool takeOverDataResponsibility) {
		if(obj!=NULL) {
			delete obj;
		}
		bool dense = true;
		size_t NNZ = 0;
		if(CheckForSparsity) {
			for(T *ptr=d,*ptr_e=d+numEl;ptr<ptr_e;++ptr) {
				NNZ += (fabs(*ptr)>1e-12) ? 1 : 0;
			}
			if(NNZ<numEl/3) {
				dense = false;
			}
		}
		if(dense) {
			if(takeOverDataResponsibility) {
				obj = new BaseClassDense<T>(d,numEl,1);
			} else {
				obj = new BaseClassDense<T>(d,numEl,0);
				/*for(size_t k=0;k<numEl;++k) {
					obj->Set(k,d[k]);
				}*/
			}
		} else {
			T* dataNew = new T[NNZ];
			size_t* ixNew = new size_t[NNZ];
			size_t* ptrIX = ixNew;
			for(T *ptr=d,*ptr_e=d+numEl,*ptrNew=dataNew;ptr<ptr_e;++ptr) {
				if(fabs(*ptr)>1e-12) {
					*ptrNew++ = *ptr;
					*ptrIX++ = ptr-d;
				}
				if(ptrIX-ixNew==int(NNZ)) {
					break;
				}
			}
			obj = new BaseClassSparse2<T>(dataNew,ixNew,NNZ,1);
			if(takeOverDataResponsibility) {
				delete [] d;
			}
		}		
	}
	ProxyObject<T> operator[](size_t ix) {
		if(obj==NULL) {
			assert(false);
		}
		return ProxyObject<T>(obj,ix);
	}
};

#endif