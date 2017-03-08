/*
 * gtf.h
 *
 *  Created on: 24.06.2015
 *      Author: kaisers
 */

#ifndef GTF_H_
#define GTF_H_

#include <fstream>
#include <string>
#include "tokenize.h"


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// First eight columns:
// 1	seqname
// 2	source
// 3	feature	e.g. "CDS", "start_codon", "stop_codon", "exon"
// 4	start 	(1-based)
// 5	end 		(inclusive)
// 6	score	(0 - 1000), No score = "."
// 7	frame	(0-2), Not a coding exon: "."
// 8	group = attributes
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

namespace gtf
{

typedef std::string::const_iterator ssci;

struct gff_element
{
public:
	size_t 		id_;
	std::string seqname_;
	std::string source_;
	std::string feature_;
	std::string start_;
	std::string end_;
	std::string score_;
	std::string strand_;
	std::string frame_;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Process first seven columns which are
	// fixed in GFF file format
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

	ssci parse_gff(ssci begin, ssci end, size_t id, char c='\t')
	{
		id_=id;
		std::string::const_iterator iter;

		// seqname
		iter = std::find(begin, end, c);
		seqname_ = std::string(begin, iter);
		if(iter == end) return end;
		begin=++iter;

		// source
		iter = std::find(begin, end, c);
		source_ = std::string(begin, iter);
		if(iter == end)	return end;
		begin=++iter;

		// feature
		iter = std::find(begin, end, c);
		feature_ = std::string(begin, iter);
		if(iter==end) return end;
		begin= ++iter;

		// start
		iter = std::find(begin, end, c);
		start_ = std::string(begin, iter);
		if(iter==end) return end;
		begin=++iter;

		// end
		iter = std::find(begin, end, c);
		end_ = std::string(begin, iter);
		if(iter==end) return end;
		begin=++iter;

		// score
		iter = std::find(begin, end, c);
		score_ = std::string(begin, iter);
		if(iter==end) return end;
		begin=++iter;

		// strand
		iter = std::find(begin, end, c);
		strand_ = std::string(begin, iter);
		if(iter==end) return end;
		begin=++iter;

		// frame
		iter=std::find(begin, end, c);
		frame_ = std::string(begin, iter);
		if(iter==end) return end;

		// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
		// Return iterator pointing at first position
		// behind last identified delimiter (c)
		// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
		return ++iter;
	}

	operator std::string() const
	{
		char delim = '\t';
		std::stringstream sst;
		sst << id_ << delim	<< "seqname"	<< delim << seqname_	<< "\n";
		sst << id_ << delim << "source"		<< delim << source_		<< "\n";
		sst << id_ << delim << "feature"	<< delim << feature_	<< "\n";
		sst << id_ << delim << "start"		<< delim << start_		<< "\n";
		sst << id_ << delim << "end"		<< delim << end_		<< "\n";
		sst << id_ << delim << "score"		<< delim << score_		<< "\n";
		sst << id_ << delim << "frame"		<< delim << frame_		<< "\n";
		return sst.str();
	}

	void process_data(void (*func)(size_t s, const char* feature, const char *value, void* o), void *o) const
	{
		func(id_, "seqname", seqname_.c_str(), o);
		func(id_, "source", source_.c_str(), o);
		func(id_, "feature", feature_.c_str(), o);
		func(id_, "start", start_.c_str(), o);
		func(id_, "end", end_.c_str(), o);
		func(id_, "score", score_.c_str(), o);
		func(id_, "frame", frame_.c_str(), o);
	}
};



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Content of the last (8th) column in GTF - file
//
// Description on UCSC site:
// - Each attribute consists of a type/value pair.
// - Attributes must end in a semi-colon, and be separated
//		from any following attribute by exactly one space.
//
// Example:
// gene_id 'ENSG00000227232'; transcript_id 'ENST00000438504';
//
// gtf_attribute parses in two levels:
// A) token_list (list of strings) separates attribute (separated by ';')
//		pairs into tokens.
// B)
//
//
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //


class gtf_attribute
{
private:

	typedef tokenize::pair_token<std::string, tokenize::wrapped_token<'"'> > wrap_pair;

	// Container
	tokenize::token_list tk_list_;
	wrap_pair wpair_;
	tokenize::token_map map_;

public:
	gtf_attribute(): tk_list_(tokenize::token_list(';')), wpair_(wrap_pair(' ')) {}
	~gtf_attribute() {}


	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Read attribute data from GTF file stream.
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	void parse(ssci begin, ssci end, size_t id)
	{
		// "gene_id 'ENSG00000227232'; transcript_id 'ENST00000438504'; ..."
		// Using delimiter ';'
		tk_list_.parse(begin, end);


		std::list<std::string>::const_iterator iter;
		for(iter=tk_list_.l_.begin(); iter != tk_list_.l_.end(); ++iter)
		{
			// "gene_id 'ENSG00000227232'"
			wpair_.parse_trim(iter->begin(), iter->end());
			// Add to unordered map
			map_.add(wpair_.string1(), id, wpair_.string2());
		}
	}


	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Data retrieval
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	size_t map_size(const std::string &map) { return map_.map_size(map); }
	operator std::string() const { return std::string(map_); }
	std::list<std::pair<size_t, std::string> > map_data() { return map_.map_data(); }

	// Generic data retrieval mechanism
	void process_data(void (*func)(size_t s, const char* feature, const char *value, void* o), void *o) const
	{ map_.process_data(func, o); }


	void process_feature(const char* feature, void (*func)(size_t s, const char* value, void *o), void *o)
	{ map_.process_feature(feature, func, o); }

};




class gtf_file
{
public:
	gtf_file(char delim = '\t'): id_(0), delim_(delim) {}
	~gtf_file() {}

private:
	size_t id_;
	char delim_;
	std::string line;
	std::string::const_iterator begin;
	std::string::const_iterator end;

	gff_element elem_;
	typedef std::list<gff_element> lge;
	lge l_;
	gtf_attribute att_;


public:
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Extract data from file.
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

	void process_lines(std::ifstream &file, char com='#')
	{
		while(getline(file, line))
		{
			if(line[0] != com)
			{
				++id_;
				begin=line.begin();
				end=line.end();
				begin = elem_.parse_gff(begin, end, id_, delim_);
				l_.push_back(elem_);
				att_.parse(begin, end, id_);
			}
		}
	}

	// Takes callback for printing out progression indicator ...
	void process_lines(std::ifstream &file, size_t nlines, void (*callback)(size_t n), char com='#')
	{
		size_t i= 0;
		while(getline(file, line))
		{
			if(!(++i % nlines))
				callback(i);

			if(line[0] != com)
			{
				++id_;
				begin=line.begin();
				end=line.end();
				begin = elem_.parse_gff(begin, end, id_);
				l_.push_back(elem_);
				att_.parse(begin, end, id_);
			}
		}
		callback(i);
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Content retrieval
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

	size_t map_size(const std::string &map) { return att_.map_size(map); }
	size_t size() const { return l_.size(); }

	std::list<std::pair<size_t, std::string> > map_data() { return att_.map_data(); }


	// Used as generic data retrieval mechanism.
	void process_data(void (*func)(size_t s, const char* feature, const char *value, void* o), void *o) const
	{
		for(lge::const_iterator gi=l_.begin(); gi!=l_.end(); ++gi)
			gi->process_data(func, o);
		att_.process_data(func, o);
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Retrieve data from first eight GTF columns:
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

	void fill_seqname(void (*func)(size_t s, const char* value, void *o), void *o)
	{
		for(lge::const_iterator gi=l_.begin(); gi!=l_.end(); ++gi)
			func(gi->id_, gi->seqname_.c_str(), o);
	}

	void fill_source(void (*func)(size_t s, const char* value, void *o), void *o)
	{
		for(lge::const_iterator gi=l_.begin(); gi!=l_.end(); ++gi)
			func(gi->id_, gi->source_.c_str(), o);
	}

	void fill_feature(void (*func)(size_t s, const char* value, void *o), void *o)
	{
		for(lge::const_iterator gi=l_.begin(); gi!=l_.end(); ++gi)
			func(gi->id_, gi->feature_.c_str(), o);
	}

	void fill_start(void (*func)(size_t s, const char* value, void *o), void *o)
	{
		for(lge::const_iterator gi=l_.begin(); gi!=l_.end(); ++gi)
			func(gi->id_, gi->start_.c_str(), o);
	}

	void fill_end(void (*func)(size_t s, const char* value, void *o), void *o)
	{
		for(lge::const_iterator gi=l_.begin(); gi!=l_.end(); ++gi)
			func(gi->id_, gi->end_.c_str(), o);
	}

	void fill_score(void (*func)(size_t s, const char* value, void *o), void *o)
	{
		for(lge::const_iterator gi=l_.begin(); gi!=l_.end(); ++gi)
			func(gi->id_, gi->score_.c_str(), o);
	}

	void fill_strand(void (*func)(size_t s, const char* value, void *o), void *o)
	{
		for(lge::const_iterator gi=l_.begin(); gi!=l_.end(); ++gi)
			func(gi->id_, gi->strand_.c_str(), o);
	}

	void fill_frame(void (*func)(size_t s, const char* value, void *o), void *o)
	{
		for(lge::const_iterator gi=l_.begin(); gi!=l_.end(); ++gi)
			func(gi->id_, gi->frame_.c_str(), o);
	}

	void proc_attr_feature(const char* feature, void (*func)(size_t s, const char* value, void *o), void *o)
	{ att_.process_feature(feature, func, o); }
};


static void print_gtf_file(size_t s, const char* feature, const char *value, void* o)
{
	std::ostream *os = (std::ostream *) o;
	*os << s << "\t" << feature << "\t" << value << "\n";
}


std::ostream & operator << (std::ostream &os, const gtf::gtf_file &gtf)
{
	gtf.process_data(print_gtf_file, &os);
	return os;
}


} // Namespace gtf


#endif /* GTF_H_ */
