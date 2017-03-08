/*
 * exptr.h
 *
 *  Created on: 30.01.2015
 *      Author: kaisers
 */

#ifndef EXPTR_HPP_
#define EXPTR_HPP_

#include<Rdefines.h>
#include<Rinternals.h>
#include<memory>
#include<vector>
#include<string>


template<typename T>
static void _finalizer(SEXP ext)
{
	if(TYPEOF(ext)==EXTPTRSXP)
	{
		std::shared_ptr<T> *p = (std::shared_ptr<T>*) R_ExternalPtrAddr(ext);
		delete p;
	}
}



template<typename T>
class extptr
{
public:

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Constructors
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	explicit extptr(): sp_(new std::shared_ptr<T>(std::make_shared<T>())){
		//Rprintf("[extptr] extptr() id: %3i\n", (*sp_)->getValue());
	}
	// Copies incoming object
	extptr(T &p): sp_(new std::shared_ptr<T>(std::make_shared<T>())) { **sp_ = p;	}


	explicit extptr(SEXP pPtr)//: sp_(new shared_ptr<T>(make_shared<T>())) // ToDo: Throw exception ??
	{
		if(TYPEOF(pPtr) != EXTPTRSXP)
			error("[extptr] No external pointer!");

		if(!R_ExternalPtrAddr(pPtr))
			error("[extptr] Received Nil pointer!");

		std::shared_ptr<T> * sp = (std::shared_ptr<T> *) (R_ExternalPtrAddr(pPtr));
		sp_ = new std::shared_ptr<T>(*sp);

	}

	virtual ~extptr()
	{
		//Rprintf("[extptr] ~exptr ptr : %p\n", sp_);
		delete sp_;
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Operator overloading
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
    extptr& operator=(const extptr& rhs)
    {
    	if(this != &rhs)
    		*sp_ = *(rhs.sp_); // Share pointer
   	    return *this ;
    }

    T& operator*() const { return *(*sp_); }
    T* operator->() const { return sp_->get(); }

	operator SEXP() const
	{
		std::shared_ptr<T> * p = new std::shared_ptr<T>;
		// Shared copy
		*p = *sp_;
		//Rprintf("[extptr] operator SEXP: %p\n", p);

	    SEXP ext = PROTECT(R_MakeExternalPtr(p, R_NilValue, R_NilValue));
	    R_RegisterCFinalizerEx(ext, _finalizer<T>, TRUE);
	    UNPROTECT(1);
	    return ext;
	}

private:
	std::shared_ptr<T> * sp_;
};



template<typename T>
SEXP to_string(const extptr<T> &p)
{
	SEXP pRes  = PROTECT(allocVector(STRSXP, 1));
	SET_STRING_ELT(pRes, 0, mkChar(std::string(*p).c_str()));
	UNPROTECT(1);
	return pRes;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// General template declaration
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
template<typename T>
class atmptr
{
	explicit atmptr(unsigned size): protected_(false), pRob(R_NilValue), p_(0) { }
	virtual ~atmptr() { if(protected_) UNPROTECT(pRob); }

    T& operator*() const { return *p_; }
    T* operator->() const { return p_; }
    T& operator [] (int i) { return p_[i]; }
    int length() const { return length(pRob); }

	operator SEXP() const { return pRob; }

private:
	bool protected_;
	SEXP pRob;
	T * p_;
};


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Template specialisation for int (LOGICAL also is INT)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

template<>
class atmptr<int>
{
public:
	explicit atmptr(unsigned size, bool set_na=false): protected_(true), pRob(R_NilValue), p_(0)
	{
		pRob = PROTECT(allocVector(INTSXP, size));
		p_ = INTEGER(pRob);
		if(set_na)
			clear();
	}

	explicit atmptr(SEXP p): protected_(false)
	{
		  if( !Rf_isInteger(p) )
		      error("[atmptr<int>] Argument must be a integer, found %s",
		            type2char(TYPEOF(p)));

		pRob=p;
		p_ = INTEGER(pRob);
	}

	virtual ~atmptr()
	{
		if(protected_)
			UNPROTECT(1);
	}

    int& operator*() const { return *p_; }
    int* operator->() const { return p_; }
    int& operator [] (int i) { return p_[i]; }
    int size() const { return length(pRob); }

    void clear()
    {
    	long i, n = length(pRob);
    	for(i=0; i < n; ++i)
    		INTEGER(pRob)[i] = NA_INTEGER;
    }
    void set_all(int i = 0) { memset(p_, i, sizeof(int) * length(pRob)); }

	operator SEXP() const { return pRob; }

private:
	bool protected_;
	SEXP pRob;
	int * p_;
};



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Template specialisation for char
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
template<>
class atmptr<char>
{
public:
	explicit atmptr(unsigned size, bool set_na=false): protected_(true), pRob(R_NilValue)
	{
		pRob = PROTECT(allocVector(STRSXP, size));
		if(set_na)
			clear();
	}

	explicit atmptr(SEXP p): protected_(false)
	{
		  if( !Rf_isString(p) )
		      error("[atmptr<char>] Argument must be a string, found %s",
		            type2char(TYPEOF(p)));

		pRob=p;
	}

	explicit atmptr(const std::vector<std::string> &v): protected_(true)
	{
		int i, n = (int) v.size();
		pRob = PROTECT(allocVector(STRSXP, n));
		for(i = 0; i < n; ++i)
			SET_STRING_ELT(pRob, i, mkChar(v[i].c_str()));
	}

	virtual ~atmptr()
	{
		if(protected_)
			UNPROTECT(1);
	}

	void get(unsigned i, std::string &s)
	{
		s.clear();
		if(i < (unsigned) LENGTH(pRob))
			s = CHAR(STRING_ELT(pRob, i));
	}

	void set(unsigned i, const std::string &s)
	{
		if(i < (unsigned) LENGTH(pRob))
			SET_STRING_ELT(pRob, i, mkChar(s.c_str()));
	}

	void set(unsigned i, const char *c)
	{
		if(i < (unsigned) LENGTH(pRob))
			SET_STRING_ELT(pRob, i, mkChar(c));
	}

	void clear()
	{
		int i, n = LENGTH(pRob);
		for(i=0; i<n; ++i)
			SET_STRING_ELT(pRob, i, NA_STRING);
	}

    int length() const { return LENGTH(pRob); }


    atmptr<int> atoi() const
	{

    	int i, n=length();
    	atmptr<int> res(n);
    	for(i=0; i < n; ++i)
    		res[n] = (int) atol(CHAR(STRING_ELT(pRob, i)));
    	return res;
	}


	operator SEXP() const { return pRob; }

private:
	bool protected_;
	SEXP pRob;
};


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Function prototypes for filling by using function pointer
// void (*set_atm_c)(unsigned i, const char* c, void *o)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
void set_atm_c(unsigned i, const char* c, void *o)
{
	atmptr<char> *p = (atmptr<char> *) o;
	p->set(i, c);
}

void set_atm_c(unsigned i, const std::string& c, void *o)
{
	atmptr<char> *p = (atmptr<char> *) o;
	p->set(i, c);

}


SEXP to_string(const std::string & s)
{
	SEXP pRes = PROTECT(allocVector(STRSXP, 1));
	SET_STRING_ELT(pRes, 0, mkChar(s.c_str()));
	UNPROTECT(1);
	return pRes;
}


#endif /* EXPTR_HPP_ */
