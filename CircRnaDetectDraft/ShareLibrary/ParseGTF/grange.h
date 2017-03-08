/*
 * grange.hpp
 *
 *  Created on: 12.02.2015
 *      Author: kaisers
 */

#ifndef GRANGE_HPP_
#define GRANGE_HPP_

#include "data_frame.h"
#include <list>

struct grange
{
	int id;
	int seqid;
	int begin;
	int end;
	// (Right) Shift which is introduced by _Unifying at _Begin position
	unsigned ub_shift;

public:
	grange(): id(++last_id), seqid(0), begin(0), end(0), ub_shift(0) {}
	grange& operator ++() { id = ++last_id; return *this; }

	grange operator ++(int)
	{
		grange tmp(*this);
		id = ++last_id;
		return tmp;
	}


private:
	static unsigned last_id;
};

unsigned grange::last_id = 0;


class grange_list
{
public:
	grange_list(): last_id(0) {}
	~grange_list() {}

	void push_back(const grange &g) { l.push_back(g); }

	operator data_frame() const
	{
		data_frame dfr( (unsigned) l.size(), 5);
		int *id = dfr.addIntColumn("id");
		int *seqid = dfr.addIntColumn("seqid");
		int *begin = dfr.addIntColumn("begin");
		int *end = dfr.addIntColumn("end");
		int *ubs = dfr.addIntColumn("ubs");

		int i;
		std::list<grange>::const_iterator iter;

		for(i = 0, iter = l.begin(); iter != l.end(); ++i, ++iter)
		{
			id[i] = iter->id;
			seqid[i] = iter->seqid;
			begin[i] = iter->begin;
			end[i] = iter->end;
			ubs[i] = iter->ub_shift;
		}
		return dfr;
	}


private:
	unsigned last_id;
	std::list<grange> l;
};




#endif /* GRANGE_HPP_ */
