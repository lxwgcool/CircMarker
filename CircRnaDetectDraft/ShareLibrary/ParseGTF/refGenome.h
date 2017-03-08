/*
 * refGenome.h
 *
 *  Created on: 25.02.2015
 *      Author: kaisers
 */

#ifndef REFGENOME_H_
#define REFGENOME_H_


// Must be included before R header because of
// collision with Rinternals length makro -> Rf_length(x)
#include <fstream>
#include "gtf.h"

#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Rdynload.h> // DllInfo

#include "data_frame.h"		// includes Rdefines.h
#include "grange.h"			// includes data_frame.h
#include "extptr.h"



extern "C" {


///////////////////////////////////////////////////////////////////////////////////////////////////
// Split Ensembl gtf Attribute column data
// and return data.frame
///////////////////////////////////////////////////////////////////////////////////////////////////
SEXP read_gtf(SEXP pParam, SEXP pProgress);

///////////////////////////////////////////////////////////////////////////////////////////////////
// Calculate exon_number from subsequent transcript and start values
// Expects ordering by transcript, seqid, start, end
///////////////////////////////////////////////////////////////////////////////////////////////////

SEXP get_exon_number(SEXP pTranscript,SEXP pSeqid, SEXP pStart, SEXP pEnd);
SEXP get_splice_juncs(SEXP pTranscript,SEXP pId,SEXP pStart,SEXP pEnd);
SEXP unify_splice_juncs(SEXP pSeqid, SEXP pLstart, SEXP pLend, SEXP pRstart, SEXP pRend,
							SEXP pId, SEXP pGeneId, SEXP pStrand, SEXP pNnmd);

///////////////////////////////////////////////////////////////////////////////////////////////////
// Removes overlaps between annotated regions
///////////////////////////////////////////////////////////////////////////////////////////////////
SEXP unify_genomic_ranges(SEXP pId, SEXP pSeqid, SEXP pBegin, SEXP pEnd);


///////////////////////////////////////////////////////////////////////////////////////////////////
// overlap: Compares list of query and reference ranges
// and reports overlaps in data.frame
///////////////////////////////////////////////////////////////////////////////////////////////////

#define OVERLAP_NO_OVER 1
#define OVERLAP_R_OVER  2
#define OVERLAP_B_OVER  3
#define OVERLAP_N_OVER  4
#define OVERLAP_L_OVER  5

SEXP overlap_ranges(SEXP qryid, SEXP qrystart, SEXP qryend, SEXP refid,SEXP refstart,SEXP refend);


///////////////////////////////////////////////////////////////////////////////////////////////////
//	Annotation of gap-sites
///////////////////////////////////////////////////////////////////////////////////////////////////

SEXP get_cum_max(SEXP pVal);
SEXP gap_overlap(SEXP pQid, SEXP pQlstart, SEXP pQlend, SEXP pQrstart, SEXP pQrend,
		SEXP pRid, SEXP pRlstart, SEXP pRlend, SEXP pRrstart, SEXP pRrend, SEXP pRmaxRend);

///////////////////////////////////////////////////////////////////////////////////////////////////
//	Routines registration
///////////////////////////////////////////////////////////////////////////////////////////////////

void R_init_refGenome(DllInfo *info);


} // extern "C"

#endif /* REFGENOME_H_ */
