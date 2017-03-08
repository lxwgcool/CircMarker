/*
 * tokenize.h
 *
 *  Created on: 15.04.2015
 *      Author: kaisers
 */

#ifndef TOKENIZE_H_
#define TOKENIZE_H_

#include <iostream>
#include <sstream>
#include <algorithm>		// transform, find_if_not
#include <list>
#include <cctype>		// isspace
#include <utility>		// pair
#include <unordered_map>

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Forward declaration: used in class stoken_list
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
namespace gtf{
	class gtf_attribute;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Global namespace for this file
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
namespace tokenize
{
typedef std::string::const_iterator ssci;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
//
//	Section: Basic functions
//				r_find_if_not, toUpper, trim, unwrap, operator<<, tokenize
//
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Returns iterator to last after skip of trailing
// characters satisfying pred (UnaryPredicate)
//
// Example:
// range [first, ..., PPPPPP, last) is transformed to
// range [first, ..., last)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

template<class InputIterator, class UnaryPredicate>
InputIterator r_find_if_not (InputIterator first, InputIterator last, UnaryPredicate pred)
{
	if(first == last)
		return first;

	--last;
	while (last!=first)
	{
		if (!pred(*last))
			return ++last;

		--last;
	}
	return first;
}

inline std::string trim(const std::string &s)
{
	std::string::const_iterator fi, bi;
	fi = std::find_if_not(s.begin(), s.end(), [](int c){return std::isspace(c);});
	bi = std::find_if_not(s.rbegin(),s.rend(),[](int c){return std::isspace(c);}).base();
	return (bi<=fi? std::string() : std::string(fi, bi));
}


inline std::string trim(ssci begin, ssci end)
{
	std::string::const_iterator fi, bi;
	fi = std::find_if_not(begin, end, [](int c){return std::isspace(c);});
	bi =    r_find_if_not(fi,    end, [](int c){return std::isspace(c);});
	return (bi<=fi? std::string() : std::string(fi, bi));
}

// Trim by moving of iterators
void ref_trim(std::string::const_iterator &b, std::string::const_iterator &e)
{
	b = std::find_if_not  (b, e, [](int c){return std::isspace(c);});
	e = r_find_if_not(b, e, [](int c){return std::isspace(c);});
}


// Returns string where leading and trailing delim - characters are removed
inline std::string unwrap(ssci begin, ssci end, char delim)
{
	while((*begin == delim) && (begin != end))
		++begin;

	// Move front before searching because 'end' points
	// *behind* last character of token
	--end;
	while((*end == delim) && (begin != end))
		--end;

	return std::string(begin, ++end);
}


//std::string toUpper(std::string s)
//{
//	std::transform(s.begin(), s.end(), s.begin(), ::toupper);
//	return s;
//}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Functional solution which fills a list
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

void tokenize(ssci begin, ssci end, std::list<std::string> &l, const char c)
{
	l.clear();
	if(begin == end)
		return;
	ssci lhs = begin, rhs;
	while(lhs != end)
	{
		rhs = find(lhs, end, c);
		l.push_back(std::string(lhs, rhs));
		if(rhs!=end)
			++rhs;
		lhs = rhs;
	}
	return;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Extracts tokens and trims strings before inserting into list
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
void trim_token(ssci begin, ssci end, std::list<std::string> &l, const char c)
{
	l.clear();
	if(begin == end)
		return;
	ssci lhs = begin, rhs;
	while(lhs != end)
	{
		rhs = find(lhs, end, c);
		l.push_back(trim(lhs, rhs));
		if(rhs!=end)
			++rhs;
		lhs = rhs;
	}
	return;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Generic function calling static function for each token
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
void tokenize(ssci begin, ssci end, void (*func)(ssci b, ssci e), const char c)
{
	ssci lhs = begin, rhs;
	while(lhs != end)
	{
		rhs = std::find(lhs, end, c);
		func(lhs, rhs);

		if(rhs != end)
			++rhs;

		lhs = rhs;
	}
	return;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Generic function calling function on some object pointer o
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
void tokenize(ssci begin, ssci end, void *o, void (*func)(ssci b, ssci e, void *o), const char c)
{
	ssci lhs = begin, rhs;
	while(lhs != end)
	{
		rhs = std::find(lhs, end, c);
		func(lhs, rhs, o);

		if(rhs != end)
			++rhs;

		lhs = rhs;
	}
	return;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
//
//	Section : Class declarations
//
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Abstract base class
// Container for tokenized strings
//
// Derived classes:
// - token-pair  (for exactly two sub-strings)
// - value_token (for type-value pairs)
// - token_list  (arbitrary number of sub-strings, uses *trim_token* function
//
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

class basic_token
{
public:
	basic_token(char delim=':'): delim_(delim) {}
	virtual ~basic_token() {}
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	std::string getDelim() const { return std::string(1, delim_); }
	void setDelim(char delim) { delim_ = delim; }
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	friend std::ostream& operator << (std::ostream&os, const basic_token &t);

protected:
	char delim_;
	virtual operator std::string() const = 0;
};

// Provides insertion operator for all derived classes
std::ostream& operator << (std::ostream &os, const basic_token &bt)
	{ return os << std::string(bt); }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
//
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
class string_token : public basic_token
{
public:
	string_token() {}
	string_token(ssci begin, ssci end): s_(begin, end) {}
	virtual ~string_token() noexcept {}

	virtual void parse(ssci begin, ssci end)
			{ s_ = std::string(begin, end); }

protected:
	virtual operator std::string() const { return s_; }

private:
	std::string s_;
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// wrapped_token:
// Contains a string variable.
// During construction from string-iterators, leading and lagging
// delim - characters (e.g. '"') are skipped.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

template<unsigned char delim>
class wrapped_token : public basic_token
{
public:
	wrapped_token() {}
	wrapped_token(ssci begin, ssci end)
	{
		while((*begin == delim) && (begin != end))
			++begin;

		// Move front before searching because 'end' points
		// *behind* last character of token
		--end;
		while((*end == delim) && (begin != end))
			--end;

		s_ = std::string(begin, ++end);
	}
	virtual ~wrapped_token() noexcept {}

private:
	std::string s_;

public:
	void clear() { s_.clear(); }
	operator std::string() const { return s_; }
};


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// pair_token
// basic_token derived class containing two template elements
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

template<typename T1, typename T2>
class pair_token: public basic_token
{
public:
	pair_token(char delim): basic_token(delim){}
	pair_token(ssci begin, ssci end, char delim): basic_token(delim)
								{ parse(begin, end); }
	virtual ~pair_token() {}

private:
	T1 first_;
	T2 second_;

protected:
	virtual operator std::string() const { return std::string(first_) + std::string(1, delim_) + std::string(second_); }


public:
	void parse(ssci begin, ssci end)
	{
		if(begin==end)
		{
			first_.clear();
			second_.clear();
			return;
		}

		// Finds first occurrence of delim
		ssci i = find(begin, end, delim_);
		if(i != end)
		{
			// Split text into first and second token
			first_ = T1(begin, i);
			second_ = T2(++i, end);
		}
		else
		{
			// Whole text goes into first token
			first_ =  T1(begin, end);
			second_.clear();
		}
		return;
	}

	// Also trims both members
	void parse_trim(ssci begin, ssci end)
	{
		if(begin==end)
		{
			first_.clear();
			second_.clear();
			return;
		}

		// Finds first occurrence of delim
		ssci i = find(begin, end, delim_);
		if(i != end)
		{
			ssci it = i;
			ref_trim(begin, it);
			// Split text into first and second token
			first_ = T1(begin, it);

			ref_trim(++i, end);
			second_ = T2(++i, end);
		}
		else
		{
			// Whole text goes into first token
			first_ =  T1(begin, end);
			second_.clear();
		}
		return;
	}

	T1 first() const { return first_; }
	T2 second() const { return second_; }

	std::string string1() { return std::string(first_); }
	std::string string2() { return std::string(second_); }
};


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// basic_token derived class containing type-value pairs
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

class value_token : public basic_token
{
public:
	value_token(ssci begin, ssci end, char delim): basic_token(delim)
	{ parse(begin, end); }
	virtual ~value_token() {}

private:
	std::string type_;
	std::string value_;

public:
	void parse(ssci begin, ssci end)
	{
		if(begin == end)
		{
			type_ = "";
			value_ = "";
			return;
		}

		// Finds first occurrence
		ssci i = find(begin, end, delim_);
		if(i != end)
		{
			type_ = std::string(begin, i);
			value_ = std::string(++i, end);
		}
		else
		{
			type_ = std::string(begin, end);
			value_ = std::string("");
		}
	}

protected:
	virtual operator std::string() const { return type_ + std::string(1, delim_) + value_; }
};


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// basic_token derived container for arbitrary number of tokens
// stored in a list.
// Segmentation of incoming string is done by *trim_token* function which
// also trims strings before insertion into list.
//
// Text example:
// "gene_id 'ENSG00000227232'; transcript_id 'ENST00000438504'; ..."
//                            ^                                ^
// See format definition:
// https://genome.ucsc.edu/FAQ/FAQformat.html#format3
// Attributes must end in a semi-colon, and be separated from
// any following attribute by exactly one space.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

class token_list: public basic_token
{
public:
	token_list(char delim): basic_token(delim) {}
	token_list(ssci begin, ssci end, char delim): basic_token(delim)
					{ tokenize(begin, end, l_, delim_); }
	virtual ~token_list() noexcept {}

private:
	std::list<std::string> l_;
	friend class gtf::gtf_attribute;

protected:
	operator std::string() const;

public:
	void clear() { l_.clear(); }
	size_t size() const { return l_.size(); }
	void parse(ssci begin, ssci end)
	{
		l_.clear();
		trim_token(begin, end, l_, delim_);
	}
};


token_list::operator std::string() const
{
	std::stringstream sst;
	if(l_.size())
	{
		std::list<std::string>::const_iterator iter = l_.begin();
		sst << *iter;
		++iter;
		for(; iter!=l_.end(); ++iter)
			sst << delim_ << *iter;
	}
	return sst.str();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Token class which allows nested string segmentation
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

/*
 * Example:
 *
 * 	string s2("texta1+texta2:textb1+textb2");
	token<token<string,'+'>, ':'> tk2(s2.begin(),s2.end());
	tk2.setDelim('#');
	cout << "[token2      ] " << tk2 << "\n";
 *
 */


template<typename T, unsigned char delim>
class token: public basic_token
{
public:
	token(): basic_token(delim) {}
	token(ssci begin, ssci end): basic_token(delim) { parse(begin, end); }

	virtual ~token()
	{
		T* t;
		while(l_.size())
		{
			t = *l_.begin();
			delete t;
			l_.pop_front();
		}
	}

private:
	std::list<T*> l_;

	// Does not clear list
	void parse(ssci begin, ssci end)
	{
		if(begin == end)
			return;

		ssci lhs = begin, rhs;
		while(lhs != end)
		{
			rhs = find(lhs, end, delim_);

			l_.push_back(new T(lhs, rhs));

			// Advance one position (move behind delim)
			if(rhs!=end)
				++rhs;

			// Do left trim
			rhs = find_if_not(rhs, end, [](int c){return std::isspace(c);});

			// Start position for next search
			lhs = rhs;
		}
	}

protected:
	virtual operator std::string() const;
};

template<typename T, unsigned char delim>
token<T, delim>::operator std::string() const
{
	std::stringstream sst;
	if(l_.size())
	{
		typename std::list<T*>::const_iterator iter = l_.begin();
		sst << *(*iter);
		++iter;
		for(; iter!=l_.end(); ++iter)
			sst << delim_ << *(*iter);
	}
	return sst.str();
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// id_token: token containing ID value
// token_map: unordered_map: Second element is id_token
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

class id_token
{
public:
	id_token() {}
	virtual ~id_token() {}

private:
	typedef std::pair<size_t, std::string> ip;
	typedef std::list<ip> lip;
	typedef std::list<ip>::const_iterator ilip;

	lip l_;
	friend class token_map;

public:
	void add(size_t id, const std::string & val)
	{
		ip i(id, val);
		l_.push_back(i);
	}

	size_t size() const { return l_.size(); }
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// token_map:
// Stores token pairs in a has-map: string + id_token
// E.g. token pair:
// transcript_id "ENST00000438504"
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
class token_map
{
public:
	token_map() {}
	virtual ~token_map() noexcept {}

	typedef std::unordered_map<std::string, id_token> umit;
	typedef umit::const_iterator imit;

private:
	umit m_;
	friend std::ostream& operator << (std::ostream&os, const token_map &t);

public:
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Insertion of data elements
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	void add(const std::string& key, size_t id, const std::string &val)
		{ m_[key.c_str()].add(id, val); }

	void add(const char *key, size_t id, const std::string &val)
		{ m_[key].add(id, val);	}

	void add(const char *key, size_t id, const char* val)
		{ m_[key].add(id, std::string(val)); }


	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Data retrieval
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	size_t size() const { return m_.size(); }

	size_t map_size(const std::string &s) { return m_[s.c_str()].size(); }

	std::list<std::pair<size_t, std::string> > map_data()
	{
		std::list<std::pair<size_t, std::string> > l;
		std::pair<size_t, std::string> p;
		umit::const_iterator iter;
		for(iter = m_.begin(); iter!=m_.end(); ++iter)
		{
			p.second=iter->first;
			p.first=iter->second.size();
			l.push_back(p);
		}
		return l;
	}


	void process_data(void (*func)(size_t s, const char* feature, const char *value, void* o), void *o) const
	{
		std::unordered_map<std::string, id_token>::const_iterator umit;
		std::list<std::pair<size_t, std::string> >::const_iterator lpit;
		for(umit = m_.begin(); umit != m_.end(); ++umit)
		{
			for(lpit=umit->second.l_.begin(); lpit!= umit->second.l_.end(); ++lpit)
				func(lpit->first, umit->first.c_str(), lpit->second.c_str(), o);
		}
	}

	void process_feature(const char* feature, void (*func)(size_t s, const char* value, void* o), void* o)
	{
		typedef std::list<std::pair<size_t, std::string> > lps;
		id_token &id = m_[feature];
		lps & l = id.l_;
		lps::const_iterator iter;
		for(iter=l.begin(); iter!=l.end(); ++iter)
			func(iter->first, iter->second.c_str(), o);
	}


	// An example showing how to extract content from object.
	operator std::string() const
	{
		char delim='\t';
		std::stringstream sst;
		std::unordered_map<std::string, id_token>::const_iterator umit;
		std::list<std::pair<size_t, std::string> >::const_iterator lpit;

		for(umit = m_.begin(); umit != m_.end(); ++umit)
		{
			for(lpit=umit->second.l_.begin(); lpit!= umit->second.l_.end(); ++lpit)
			{
				sst << lpit->first << delim << umit->first << delim << lpit->second << "\n";
			}
		}
		return sst.str();
	}
};


std::ostream & operator << (std::ostream &os, const token_map &t)
{
	std::unordered_map<std::string, id_token>::const_iterator iter;

	for(iter = t.m_.begin(); iter != t.m_.end(); ++iter)
	{
			os << iter->first << ": " << iter->second.size() << "\n";
	}
	return os;
}

} // End namespace token

#endif /* TOKENIZE_H_ */
