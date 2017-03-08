/*
 * refGenome.cpp
 *
 *  Created on: 25.02.2015
 *      Author: kaisers
 */

#ifndef REFGENOME_CPP_
#define REFGENOME_CPP_


#include "refGenome.h"


extern "C" {

static const int buf_size=2048; // buffer size for printing ints into chars

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Calculate exon_number from subsequent transcript and start values
// Expects ordering by transcript, seqid, start, end
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

SEXP get_splice_juncs(SEXP pTranscript,SEXP pId, SEXP pStart, SEXP pEnd)
{
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Expects that pTranscript is ordered
	// so that consecutive equal pTranscript values
	// represent junctions
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Check incoming args

	if( !Rf_isInteger(pTranscript) )
		error("first argument must be a integer, found %s",
				type2char(TYPEOF(pTranscript)));

	if( !Rf_isInteger(pId) )
		error("second argument must be a integer, found %s",
				type2char(TYPEOF(pId)));

	if( !Rf_isInteger(pStart) )
		error("third argument must be a integer, found %s",
				type2char(TYPEOF(pStart)));

	if( !Rf_isInteger(pEnd) )
		error("fourth argument must be a integer, found %s",
				type2char(TYPEOF(pEnd)));

	unsigned inRow = (unsigned) LENGTH(pTranscript);
	if(((unsigned) LENGTH(pId) != inRow) | ((unsigned)LENGTH(pStart)!=inRow) | ((unsigned)LENGTH(pEnd)!=inRow))
		error("[get_splice_juncs] All arguments must have same length!");


	int *tr=INTEGER(pTranscript), *id=INTEGER(pId);
	int *start=INTEGER(pStart), *end=INTEGER(pEnd);

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Count number of splice junctions
	// = Number of row pairs where tr[i]==tr[i+1]
	unsigned nJunc=0, i,j;
	for(i=1, j=0; i < inRow; ++i, ++j)
	{
		if(tr[i]==tr[j])
			++nJunc;
	}

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Create output vectors
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	SEXP pLexid;
	PROTECT(pLexid=allocVector(INTSXP,nJunc));
	SEXP pRexid;
	PROTECT(pRexid=allocVector(INTSXP,nJunc));
	SEXP pLstart;
	PROTECT(pLstart=allocVector(INTSXP,nJunc));
	SEXP pLend;
	PROTECT(pLend=allocVector(INTSXP,nJunc));
	SEXP pRstart;
	PROTECT(pRstart=allocVector(INTSXP,nJunc));
	SEXP pRend;
	PROTECT(pRend=allocVector(INTSXP,nJunc));

	int *lexid=INTEGER(pLexid);
	int *rexid=INTEGER(pRexid);
	int *lstart=INTEGER(pLstart);
	int *lend=INTEGER(pLend);
	int *rstart=INTEGER(pRstart);
	int *rend=INTEGER(pRend);

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Calculate splice pairs
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	unsigned iJunc=0;
	for(i=1,j=0;i<inRow;++i,++j)
	{
		if(iJunc>nJunc)
			error("[get_splice_juncs] iJunc error: i=%i\tnJunc=%i\n",i,nJunc);

		if(tr[i]==tr[j])
		{
			lexid[iJunc]=id[j];
			lstart[iJunc]=start[j];
			lend[iJunc]=end[j];

			rexid[iJunc]=id[i];
			rstart[iJunc]=start[i];
			rend[iJunc]=end[i];
			++iJunc;
		}
	}

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Assemble result
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,6));
	SET_VECTOR_ELT(dflist,0,pLexid);
	SET_VECTOR_ELT(dflist,1,pRexid);
	SET_VECTOR_ELT(dflist,2,pLstart);
	SET_VECTOR_ELT(dflist,3,pLend);
	SET_VECTOR_ELT(dflist,4,pRstart);
	SET_VECTOR_ELT(dflist,5,pRend);

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Column Names
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,6));
	SET_STRING_ELT(col_names,0,mkChar("lexid"));
	SET_STRING_ELT(col_names,1,mkChar("rexid"));
	SET_STRING_ELT(col_names,2,mkChar("lstart"));
	SET_STRING_ELT(col_names,3,mkChar("lend"));
	SET_STRING_ELT(col_names,4,mkChar("rstart"));
	SET_STRING_ELT(col_names,5,mkChar("rend"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Row Names
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	SEXP row_names;
    PROTECT(row_names=allocVector(STRSXP,nJunc));

	char *buf=(char*) calloc(buf_size,sizeof(char));
    for(i=0;i<nJunc;++i)
    {
    	sprintf(buf,"%i",i+1);
    	SET_STRING_ELT(row_names,i,mkChar(buf));
    }
    free(buf);
    setAttrib(dflist,R_RowNamesSymbol,row_names);

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));

	UNPROTECT(9);
	return dflist;
}

SEXP unify_splice_juncs(SEXP pSeqid, SEXP pLstart, SEXP pLend, SEXP pRstart, SEXP pRend, SEXP pId, SEXP pGeneId, SEXP pStrand, SEXP pNnmd)
{
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Calculates unique junction coordinate values (uJunc)
	// from a sorted junction list
	//
	// Expects that pSeqid, pLend and pRstart are ordered
	// so that identical junctions appear as
	// consecutive equal values
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //


	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Check incoming arguments
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	if(TYPEOF(pSeqid) != INTSXP)
		error("[unify_splice_juncs] pSeqid must be INT!");
	if(TYPEOF(pLstart) != INTSXP)
		error("[unify_splice_juncs] pLstart must be INT!");
	if(TYPEOF(pLend) != INTSXP)
		error("[unify_splice_juncs] pLend must be INT!");
	if(TYPEOF(pRstart) != INTSXP)
		error("[unify_splice_juncs] pRstart must be INT!");
	if(TYPEOF(pRend) != INTSXP)
		error("[unify_splice_juncs] pRend must be INT!");
	if(TYPEOF(pId) != INTSXP)
		error("[unify_splice_juncs] pId must be INT!");
	if(TYPEOF(pGeneId) != INTSXP)
		error("[unify_splice_juncs] pGeneId must be INT!");
	if(TYPEOF(pStrand) != INTSXP)
		error("[unify_splice_juncs] pStrand must be INT!");
	if(TYPEOF(pNnmd) != INTSXP)
		error("[unify_splice_juncs] pNnmd must be INT!");


	unsigned i,j,k, nJunc, n = LENGTH(pId), nSites;

	int 	*seqid=INTEGER(pSeqid),
			*lstart=INTEGER(pLstart), *lend=INTEGER(pLend),
			*rstart=INTEGER(pRstart), *rend=INTEGER(pRend),
			*id=INTEGER(pId),
			*gid=INTEGER(pGeneId),    *strand=INTEGER(pStrand);
	int		*nnmd = INTEGER(pNnmd);

	int usq,ule,urs; // unified splice coordinates

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// First row identifies first nJunc site
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	usq = seqid[0];
	ule = lend[0];
	urs = rstart[0];
	nJunc = 1;
	i = 0;
	//Rprintf("[unify_splice_juncs] Start junc %i at i= %i.\n",nJunc,i);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Check all subsequent rows for position equality
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	for(i = 1; i < n; ++i)
	{
		if((seqid[i] != usq) | (lend[i] != ule) | (rstart[i] != urs))
		{
			++nJunc;
			//Rprintf("[unify_splice_juncs] Start junc %i at i= %i.\n",nJunc,i);
			usq = seqid[i];
			ule = lend[i];
			urs = rstart[i];
		}
	}
	//Rprintf("[unify_splice_juncs] Found %i juncs.\n",nJunc);

	if(nJunc == 0)
		return R_NilValue;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Create output vectors
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	SEXP puId;
	PROTECT(puId = allocVector(INTSXP, nJunc));

	SEXP puSeq;
	PROTECT(puSeq = allocVector(INTSXP, nJunc));

	SEXP puLstart;
	PROTECT(puLstart = allocVector(INTSXP, nJunc));

	SEXP puLend;
	PROTECT(puLend = allocVector(INTSXP, nJunc));

	SEXP puRstart;
	PROTECT(puRstart = allocVector(INTSXP, nJunc));

	SEXP puRend;
	PROTECT(puRend = allocVector(INTSXP, nJunc));

	SEXP puNsites; // Number of junction sites per uJunc
	PROTECT(puNsites = allocVector(INTSXP, nJunc));

	SEXP puGene; // Gene-id associated with exon table
	PROTECT(puGene = allocVector(INTSXP, nJunc));

	SEXP puStrand;
	PROTECT(puStrand = allocVector(INTSXP, nJunc));

	SEXP pFexid; // First exon id (points into exon table)
	PROTECT(pFexid = allocVector(INTSXP, nJunc));

	// Number of transcripts which are not of transcript biotype "nonsense_mediated_decay"
	// -> Count not NMD = cnnmd
	SEXP pCnNmd;
	PROTECT(pCnNmd = allocVector(INTSXP, nJunc));


	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Calculate values for sites
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	int     *uid = INTEGER(puId),         *useq = INTEGER(puSeq),
			*ulstart = INTEGER(puLstart), *ulend = INTEGER(puLend),
			*urstart = INTEGER(puRstart), *urend = INTEGER(puRend),
			*uNsites = INTEGER(puNsites), *uId = INTEGER(pFexid),
			*uGene = INTEGER(puGene),     *uStrand = INTEGER(puStrand);
	int * cnnmd = INTEGER(pCnNmd);

	int min_lstart,max_rend;
	const unsigned nGenes = 10;
	unsigned geneId[nGenes];  // geneId
	unsigned geneCt[nGenes];  // gene-count
	unsigned geneStr[nGenes]; // gene-strand
	unsigned mxGct,mxGid,mxStr;    // maxGeneCount, maxGeneId, maxStrand



	j=0; // write index of actual uJunc
	i=0; // read  index of actual junc

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Start reading in row 0
	// Row 0 identifies first uJunc
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	uid[j] = j + 1;
	useq[j] = seqid[i];
	ulend[j] = lend[i];
	urstart[j] = rstart[i];
	uId[j] = id[i];
	min_lstart = lstart[i];
	max_rend = rend[i];
	nSites = 1;
	cnnmd[j] = nnmd[i];
	//Rprintf("[unify_splice_juncs] Start junc %i at i= %i.\n",j,i);

	geneId[0] = gid[i];
	geneStr[0] = strand[i];
	geneCt[0] = 1;
	for(k = 1; k < nGenes; ++k)
		geneId[k] = 0;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Check all subsequent rows for position equality
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	for(i = 1; i < n; ++i)
	{
		if( (seqid[i] != useq[j]) | (lend[i] != ulend[j]) | (rstart[i] != urstart[j]) )
		{
			// Row i contains new uJunc-site

			// - - - - - - - - - - - - - - - - - - - - - //
			// Complete values for last uJunc-site
			// - - - - - - - - - - - - - - - - - - - - - //
			uNsites[j] = nSites;
			ulstart[j] = min_lstart;
			urend[j] = max_rend;

			// - - - - - - - - - - - - - - - - - - - - - //
			// Get some gene-id with "maximal" gene-id count
			// - - - - - - - - - - - - - - - - - - - - - //
			mxGct = 0;
			mxGid = 0;
			mxStr = 0;
			for(k = 0; k < nGenes; ++k)
			{
				if(geneId[k] == 0)
					break;

				if(mxGct < geneCt[k])
				{
					mxGct = geneCt[k];
					mxGid = geneId[k];
					mxStr = geneStr[k];
				}
			}
			uGene[j] = mxGid;
			uStrand[j] = mxStr;

			//Rprintf("[unify_splice_juncs] Start junc %i at i= %i.\n",j,i);

			// - - - - - - - - - - - - - - - - - - - - - //
			// Set values for next uJunc-site
			// - - - - - - - - - - - - - - - - - - - - - //
			++j;
			if(j >= nJunc)
				error("[unify_splice_juncs] Write index exceeding nJunc limit!");

			uid[j] = j + 1;
			useq[j] = seqid[i];
			ulend[j] = lend[i];
			urstart[j] = rstart[i];
			uId[j] = id[i];
			min_lstart = lstart[i];
			max_rend = rend[i];
			nSites = 1;
			cnnmd[j] = nnmd[i];

			geneId[0] = gid[i];
			geneStr[0] = strand[i];
			geneCt[0] = 1;

			for(k = 1; k < nGenes; ++k)
			{
				geneId[k] = 0;
				geneCt[k] = 0;
			}

		}
		else
		{
			// - - - - - - - - - - - - - - - - - - - - - //
			// Row i is part of actual uJunc-site
			// - - - - - - - - - - - - - - - - - - - - - //
			++nSites;

			// lstart and rend
			if(lstart[i] < min_lstart)
				min_lstart = lstart[i];

			if(rend[i] > max_rend)
				max_rend = rend[i];

			// Count geneId's
			for(k = 0; k < nGenes; ++k)
			{
				if((unsigned) gid[i] == geneId[k])
				{
					++(geneCt[k]);
					break;
				}
				if(geneId[k] == 0)
				{
					geneId[k] = gid[i];
					geneStr[k] = strand[i];
					geneCt[k] = 1;
					break;
				}
			}

			// - - - - - - - - - - - - - - - - - - - - - //
			// Is row i part of non-NMD transcript?
			// - - - - - - - - - - - - - - - - - - - - - //
			cnnmd[j] += nnmd[i];
		}
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Complete values for last uJunc-site
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	uNsites[j] = nSites;
	ulstart[j] = min_lstart;
	urend[j] = max_rend;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Get some gene-id with "maximal" gene-id count
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	mxGct = 0;
	mxGid = 0;
	mxStr = 0;
	for(k = 0; k < nGenes; ++k)
	{
		if(geneId[k]==0)
			break;

		if(mxGct < geneCt[k])
		{
			mxGct = geneCt[k];
			mxGid = geneId[k];
			mxStr = geneStr[k];
		}
	}
	uGene[j] = mxGid;
	uStrand[j] = mxStr;


	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Assemble output
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	SEXP dflist;
	PROTECT(dflist = allocVector(VECSXP,11));

	SET_VECTOR_ELT(dflist, 0, puId);
	SET_VECTOR_ELT(dflist, 1, puSeq);
	SET_VECTOR_ELT(dflist, 2, puLstart);
	SET_VECTOR_ELT(dflist, 3, puLend);
	SET_VECTOR_ELT(dflist, 4, puRstart);
	SET_VECTOR_ELT(dflist, 5, puRend);
	SET_VECTOR_ELT(dflist, 6, puNsites);
	SET_VECTOR_ELT(dflist, 7, puGene);
	SET_VECTOR_ELT(dflist, 8, puStrand);
	SET_VECTOR_ELT(dflist, 9, pFexid);
	SET_VECTOR_ELT(dflist, 10, pCnNmd);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Column Names
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,11));
	SET_STRING_ELT(col_names, 0, mkChar("id"));
	SET_STRING_ELT(col_names, 1, mkChar("seqid"));
	SET_STRING_ELT(col_names, 2, mkChar("lstart"));
	SET_STRING_ELT(col_names, 3, mkChar("lend"));
	SET_STRING_ELT(col_names, 4, mkChar("rstart"));
	SET_STRING_ELT(col_names, 5, mkChar("rend"));
	SET_STRING_ELT(col_names, 6, mkChar("nSites"));
	SET_STRING_ELT(col_names, 7, mkChar("gene_id"));
	SET_STRING_ELT(col_names, 8, mkChar("strand"));
	SET_STRING_ELT(col_names, 9, mkChar("fexid"));
	SET_STRING_ELT(col_names, 10, mkChar("cnnmd"));
	setAttrib(dflist, R_NamesSymbol, col_names);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Row Names
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	SEXP row_names;
    PROTECT(row_names = allocVector(STRSXP, nJunc));

	char *buf=(char*) calloc(buf_size,sizeof(char));
    for(i = 0; i < nJunc; ++i)
    {
    	sprintf(buf,"%i", i + 1);
    	SET_STRING_ELT(row_names, i, mkChar(buf));
    }
    free(buf);
    setAttrib(dflist,R_RowNamesSymbol, row_names);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	setAttrib(dflist, R_ClassSymbol, mkString("data.frame"));

	UNPROTECT(14);
	return dflist;
}


SEXP get_exon_number(SEXP pTranscript,SEXP pSeqid, SEXP pStart, SEXP pEnd)
{
	if(TYPEOF(pTranscript)!=INTSXP)
		error("[get_exon_number] pTranscript must be INT!");
	if(TYPEOF(pSeqid)!=INTSXP)
		error("[get_exon_number] pSeqid must be INT!");
	if(TYPEOF(pStart)!=INTSXP)
		error("[get_exon_number] pStart must be INT!");
	if(TYPEOF(pEnd)!=INTSXP)
		error("[get_exon_number] pEnd must be INT!");

	int n=LENGTH(pTranscript);
	if(LENGTH(pSeqid)!=n || LENGTH(pSeqid)!=n || LENGTH(pStart)!=n || LENGTH(pEnd)!=n)
		error("[get_exon_number] All args must have same length!");

	SEXP res;
	PROTECT(res=allocVector(INTSXP,n));
	int i,j, exon_number=1, nSeqMm=0, nStartMm=0;

	INTEGER(res)[0]=exon_number;

	for(i=1,j=0;i<n;++i,++j)
	{
		if(INTEGER(pTranscript)[j]==INTEGER(pTranscript)[i])
		{
			// same transcript
			++exon_number;

			// Security checks
			if(INTEGER(pSeqid)[j]!=INTEGER(pSeqid)[i])
				++nSeqMm;
			if(INTEGER(pEnd)[j]>=INTEGER(pStart)[i])
				++nStartMm;
		}
		else
			exon_number=1; // new transcript

		INTEGER(res)[i]=exon_number;
	}

	if(nSeqMm>0)
		Rprintf("[get_exon_number] Found %i sequence mismatches!\n",nSeqMm);
	if(nStartMm>0)
		Rprintf("[get_exon_number] Found %i start-end mismatches!\n",nStartMm);

	UNPROTECT(1);
	return res;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Removes overlaps between annotated regions
///////////////////////////////////////////////////////////////////////////////////////////////////
SEXP unify_genomic_ranges(SEXP pId, SEXP pSeqid, SEXP pBegin, SEXP pEnd)
{
	if(TYPEOF(pId) != INTSXP)
		error("pId must be Integer!");

	if(TYPEOF(pSeqid) != INTSXP)
		error("pSeqid must be Integer!");

	if(TYPEOF(pBegin) != INTSXP)
		error("pBegin must be Integer!");

	if(TYPEOF(pEnd) != INTSXP)
		error("pEnd must be Integer!");

	int i, n = length(pId);

	if(length(pSeqid) != n)
		error("pId and pSeqid must have equal length!");

	if(length(pBegin) != n)
		error("pId and pBegin must have equal length!");

	if(length(pEnd) != n)
		error("pId and pEnd must have equal length!");

	int *id = INTEGER(pId);
	int *seqid = INTEGER(pSeqid);
	int *begin = INTEGER(pBegin);
	int *end = INTEGER(pEnd);


	grange_list l;
	grange last;
	last.seqid = 0;

	i=0;
	last.id = id[i];
	last.seqid = seqid[i];
	last.begin = begin[i];
	last.end = end[i];
	l.push_back(last);

	for(++i; i < n; ++i)
	{
		// New seqid
		if(seqid[i] > last.seqid)
		{
			last.id = id[i];
			last.seqid = seqid[i];
			last.begin = begin[i];
			last.end = end[i];
			last.ub_shift = 0;
			l.push_back(last);
		}else{
			// No overlap
			if(begin[i] > last.end)
			{
				last.id = id[i];
				last.begin = begin[i];
				last.end = end[i];
				last.ub_shift = 0;
				l.push_back(last);

			}else{
				// Partial overlap:
				// Introduces shift of begin position
				if(end[i] > last.end)
				{
					last.id = id[i];
					last.ub_shift = last.begin;
					last.begin = last.end + 1;
					last.ub_shift = last.begin - last.ub_shift;
					last.end = end[i];
					l.push_back(last);
				}
			}
		}
	}

	return data_frame(l);
}




// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// overlap: Compares list of query and reference ranges
// and reports overlaps in data.frame
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

SEXP overlap_ranges(SEXP qryid, SEXP qrystart, SEXP qryend, SEXP refid,SEXP refstart,SEXP refend)
{
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// expects qryend >= qrystart, refend >= refstart
	// expects qrystart and refstart ascending sorted
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Check type of incoming args
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	if(TYPEOF(qryid)!=INTSXP)
		error("[overlap_ranges] qryid is no INT!\n");
	if(TYPEOF(qrystart)!=INTSXP)
		error("[overlap_ranges] qrystart is no INT!\n");
	if(TYPEOF(qryend)!=INTSXP)
		error("[overlap_ranges] qryend is no INT!\n");
	if(TYPEOF(refid)!=INTSXP)
		error("[overlap_ranges] refid is no INT!\n");
	if(TYPEOF(refstart)!=INTSXP)
		error("[overlap_ranges] refstart is no INT!\n");
	if(TYPEOF(refend)!=INTSXP)
		error("[overlap_ranges] refend is no INT!\n");

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Check size of incoming args
	// Size of qry determines size of output parameters
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	unsigned nRows=LENGTH(qryid);
	if(nRows==0)
		error("[overlap_ranges] qryid had length zero!");
	if( (unsigned) LENGTH(qrystart)!=nRows)
		error("[overlap_ranges] length(qrystart)!=length(qryid)!");
	if( (unsigned) LENGTH(qryend)!=nRows)
		error("[overlap_ranges] length(qryend)!=length(qryid)!");

	unsigned nRef=LENGTH(refid);
	if(nRef==0)
		error("[overlap_ranges] refid has length zero!");
	if((unsigned) LENGTH(refstart)!=nRef)
		error("[overlap_ranges] length(refstart)!=length(refid)!");
	if((unsigned) LENGTH(refend)!=nRef)
		error("[overlap_ranges] length(refstart)!=length(refend)!");

	// read args
	unsigned *qid   =(unsigned*)INTEGER(qryid);
	unsigned *qstart=(unsigned*)INTEGER(qrystart);
	unsigned *qend  =(unsigned*)INTEGER(qryend);
	unsigned *rid   =(unsigned*)INTEGER(refid);
	unsigned *rstart=(unsigned*)INTEGER(refstart);
	unsigned *rend  =(unsigned*)INTEGER(refend);

	unsigned nProtected=0;
	unsigned qidx=0, ridx=0; 	// query and ref indices + max idices

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Column 0: overlap code
	SEXP vov;
	PROTECT(vov=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 1: difference to annotated position
	SEXP vldiff;
	PROTECT(vldiff=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 2: difference to annotated position
	SEXP vrdiff;
	PROTECT(vrdiff=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 3: id of query item
	SEXP vqid;
	PROTECT(vqid=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 4: id of ref item
	SEXP vrid;
	PROTECT(vrid=allocVector(INTSXP,nRows));
	++nProtected;
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //


	while((qidx<nRows) & (ridx<nRef))
	{
		// qry misses right
		if(qstart[qidx] > rend[ridx])
		{
			++ridx;
			continue;
		}

		if(qend[qidx] >rend[ridx])
		{
			if(qstart[qidx] >= rstart[ridx])
			{
				// qry overhang right: R_OVER
				INTEGER(vqid)  [qidx]=qid[qidx];
				INTEGER(vrid)  [qidx]=rid[ridx];
				INTEGER(vov)   [qidx]=OVERLAP_R_OVER;
				INTEGER(vldiff)[qidx]=qstart[qidx]-rstart[ridx];
				INTEGER(vrdiff)[qidx]=qend[qidx]-rend[ridx];
				++qidx;
				continue;
			}
			else
			{
				// qry overhang on both sides: B_OVER
				INTEGER(vqid)  [qidx]=qid[qidx];
				INTEGER(vrid)  [qidx]=rid[ridx];
				INTEGER(vov)   [qidx]=OVERLAP_B_OVER;
				INTEGER(vldiff)[qidx]=rstart[ridx]-qstart[qidx];
				INTEGER(vrdiff)[qidx]=qend[qidx]-rend[ridx];
				++qidx;
				continue;
			}
		}
		else if(qend[qidx] >= rstart[ridx])
		{
			if(qstart[qidx] >= rstart[ridx])
			{
				// no qry-overhang: N_OVER
				INTEGER(vqid)  [qidx]=qid[qidx];
				INTEGER(vrid)  [qidx]=rid[ridx];
				INTEGER(vov)   [qidx]=OVERLAP_N_OVER;
				INTEGER(vldiff)[qidx]=qstart[qidx]-rstart[ridx];
				INTEGER(vrdiff)[qidx]=rend[ridx]-qend[qidx];
				++qidx;
				continue;
			}
			else
			{
				// qry overhang left: L_OVER
				INTEGER(vqid)  [qidx]=qid[qidx];
				INTEGER(vrid)  [qidx]=rid[ridx];
				INTEGER(vov)   [qidx]=OVERLAP_L_OVER;
				INTEGER(vldiff)[qidx]=rstart[ridx]-qstart[qidx];
				INTEGER(vrdiff)[qidx]=rend[ridx]-qend[qidx];
				++qidx;
				continue;
			}
		}
		else
		{
			// No overlap
			INTEGER(vqid) [qidx]=qid[qidx];
			INTEGER(vrid) [qidx]=0;
			INTEGER(vov)  [qidx]=OVERLAP_NO_OVER;
			// rdiff gives distance to next ref on right side
			INTEGER(vrdiff)[qidx]=rstart[ridx]-qend[qidx];

			// ldiff gives distance to next ref on left side (or 0 if not exists)
			if(ridx>0)
				INTEGER(vldiff)[qidx]=qstart[qidx]-rend[ridx-1];
			else
				INTEGER(vldiff)[qidx]=0;

			++qidx;
			continue;
		}
	}

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Process remaining qry rows when rightmost ref ranges are passed
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	if((qidx<nRows) & (qstart[qidx] > rend[nRef-1]))
	{
		unsigned refend=rend[nRef-1];
		while(qidx<nRows)
		{
			INTEGER(vqid)  [qidx]=qid[qidx];
			INTEGER(vrid)  [qidx]=0;
			INTEGER(vov)   [qidx]=OVERLAP_NO_OVER;
			INTEGER(vrdiff)[qidx]=0;
			INTEGER(vldiff)[qidx]=qstart[qidx]-refend;
			++qidx;
		}
	}

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Convert vov to factor
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	SEXP levs;
	int nLevels=5;
	PROTECT(levs=allocVector(STRSXP,nLevels));
	++nProtected;

	SET_STRING_ELT(levs,0,mkChar("no"));
	SET_STRING_ELT(levs,1,mkChar("r"));
	SET_STRING_ELT(levs,2,mkChar("b"));
	SET_STRING_ELT(levs,3,mkChar("n"));
	SET_STRING_ELT(levs,4,mkChar("l"));
	setAttrib(vov,R_LevelsSymbol,levs);

	SEXP csymb;
	PROTECT(csymb=mkString("factor"));
	++nProtected;
	setAttrib(vov,R_ClassSymbol,csymb);

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Create data.frame
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	unsigned nCols=5;
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;

	SET_VECTOR_ELT(dflist,0,vov);
	SET_VECTOR_ELT(dflist,1,vldiff);
	SET_VECTOR_ELT(dflist,2,vrdiff);
	SET_VECTOR_ELT(dflist,3,vqid);
	SET_VECTOR_ELT(dflist,4,vrid);

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Column Names
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(col_names,0,mkChar("overlap"));
	SET_STRING_ELT(col_names,1,mkChar("leftDiff"));
	SET_STRING_ELT(col_names,2,mkChar("rightDiff"));
	SET_STRING_ELT(col_names,3,mkChar("queryid"));
	SET_STRING_ELT(col_names,4,mkChar("refid"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Create row names for data.frame
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	SEXP row_names;
    PROTECT(row_names=allocVector(STRSXP,nRows));
    ++nProtected;

	int buf_size=1024;
	char *buf=(char*) calloc(buf_size,sizeof(char));
    for(qidx=0;qidx<nRows;++qidx)
    {
    	sprintf(buf,"%i",qidx);
    	SET_STRING_ELT(row_names,qidx,mkChar(buf));
    }
    free(buf);
    setAttrib(dflist,R_RowNamesSymbol,row_names);

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
    // Make list to data.frame
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
	UNPROTECT(nProtected);

	return dflist;
}

/*
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
//                                                                                                  #
// Split gtf Attribute column data                                                                  #
// and return list with two data.frames                                                             #
//                                                                                                  #
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
// A: The Brent Lab (Washington University, St.Louis)                                               #
// http://mblab.wustl.edu/GTF2.html                                                                 #
// [attributes] All four features have the same two mandatory attributes at the end of the record:  #
//                                                                                                  #
// gene_id value;       A globally unique identifier for the genomic source of the transcript       #
// transcript_id value; A globally unique identifier for the predicted transcript.                  #
//                                                                                                  #
// These attributes are designed for handling multiple transcripts from the same genomic region.    #
// Any other attributes or comments must appear after these two and will be ignored.                #
//                                                                                                  #
// Attributes must end in a semicolon which must then be separated from the start of any subsequent #
// attribute by exactly one space character (NOT a tab character). Textual attributes *should* be   #
// surrounded by doublequotes.                                                                      #
//                                                                                                  #
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
// B: Wellcome Trust Sanger Institute                                                               #
// http://www.sanger.ac.uk/resources/software/gff/spec.html                                         #
//                                                                                                  #
// Free text values *must* be quoted with double quotes.                                            #
// Note: all non-printing characters in such free text value strings (e.g. newlines, tabs, control  #
// characters, etc) must be explicitly represented by their C (UNIX) style backslash-escaped        #
// representation (e.g. newlines as '\n', tabs as '\t').                                            #
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#


SEXP split_gtf_attr(SEXP id_vec,SEXP attr_vec)
{
	// Check type of incoming args
	if(TYPEOF(id_vec)!=INTSXP)
		error("[split_gtf_attr] id_vec is no INT!\n");
	if(TYPEOF(attr_vec)!=STRSXP)
		error("[split_gtf_attr] attr_vec is no STR: %i!\n",TYPEOF(attr_vec));

	unsigned long int i,n;

	n=LENGTH(id_vec);
	if( (unsigned) LENGTH(attr_vec) != n)
		error("[split_gtf_attr] id_vec and attr_vec must have same length!\n");

	unsigned nProtected=0;
	int *id=INTEGER_POINTER(id_vec);

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Column 1: id vector
	SEXP idvec;
	PROTECT(idvec=allocVector(INTSXP,n));
	++nProtected;

	// Column 2: gene_id
	SEXP geneIdVec;
	PROTECT(geneIdVec=allocVector(STRSXP,n));
	++nProtected;

	// Column 3: transcript_id
	SEXP transIdVec;
	PROTECT(transIdVec=allocVector(STRSXP,n));
	++nProtected;
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //


	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// flag characters
	const char space=' ', delim=';', quote='"', zero='\0';
	// iterator positions
	const char *token_first,*token_second,*iter;
	// string length
	unsigned long first_len,second_len;
	// attribute index
	unsigned attr_index;
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

	ptr_pair_list *l=ptr_pair_list_init();
	for(i = 0; i < n; ++i)
	{
		INTEGER(idvec)[i]=id[i];
		iter=CHAR(STRING_ELT(attr_vec,i));
		attr_index=1;

		// Skip leading spaces
		while((*iter==space))
			++iter;

		while(*iter != zero)
		{
			// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
			// First token
			// Take start position
			// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
			token_first=iter;
			// Rprintf("[split_gtf_attr] token_first: '%s'\n",token_first);

			if((*iter != 'g') && (attr_index == 1))
				error("[split_gtf_attr] First item must be 'gene_id': '%s'!", iter);

			// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
			// Ensembl has introduced gene records with version 76.
			// Their second attribute is 'gene_name' (unlike 'transcript_id')
			// These gene_name records are silently accepted.
			// read.gtf later splits the table when gene records are present.
			// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
			if((*iter != 't') && (attr_index == 2))
			{
				if( (*iter != 't') && (*iter != 'g'))
				{
					Rprintf("[split_gtf_attr] In line %i: '%s'\n", i + 1, CHAR(STRING_ELT(attr_vec,i)));
					error("[split_gtf_attr] Second item must be 'transcript_id' or 'gene_name': '%s'!", iter);
				}
			}

			// proceed until space and take length
			while((*iter!=space) && (*iter!=zero))
				++iter;

			if(*iter == zero)
				error("[split_gtf_attr] Found end of string in first token in line %lu: '%s'!\n",i + 1, iter);


			if(iter == token_first)
			{
				Rprintf("[split_gtf_attr] In line %i: '%s'\n", i + 1, CHAR(STRING_ELT(attr_vec,i)));
				error("[split_gtf_attr] First token ist empty: '%s'!\n",iter);
			}

			first_len = iter - token_first;

			// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
			// second token:
			// skip spaces and quotes, then take start position
			// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

			while((*iter == space) || (*iter == quote))
				++iter;

			token_second=iter;
			// Rprintf("[split_gtf_attr] token_second: '%s'\n",token_second);

			// proceed until space or quote
			while((*iter!=space) && (*iter!=quote) && (*iter!=delim) && (*iter!=zero))
				++iter;
			second_len=iter-token_second;

			// second token may be closed by quote
			if(*iter==quote)
				++iter;

			// second token must be terminated by delim
			if(*iter!=delim)
			{
				Rprintf("[split_gtf_attr] In line %i: '%s'\n", i + 1, CHAR(STRING_ELT(attr_vec,i)));
				error("[split_gtf_attr] Second token must end on ';': '%s'!",token_second);
			}
			++iter;

			// There may be terminating spaces
			while(*iter==space)
				++iter;

			// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
			//  Process collected pointer positions
			//  First token must be gene_id
			// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
			if(attr_index==1)
				SET_STRING_ELT(geneIdVec,i,mkCharLen(token_second,second_len));
			// Second token must be transcript_id
			else if(attr_index==2)
				SET_STRING_ELT(transIdVec,i,mkCharLen(token_second,second_len));
			// Probably some more tokens
			else
				ptr_pair_list_push_back(l,token_first,first_len,token_second,second_len,id[i]);
			++attr_index;
		}
	}

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	//  Create output data.frames
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	unsigned nCols, nRows;
	char buf[20];

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Merge id,gene_id and transcript_id into idflist data.frame
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	nCols = 3;
	nRows = (unsigned) n;
	SEXP idflist;
	PROTECT(idflist = allocVector(VECSXP,nCols));
	++nProtected;

	SET_VECTOR_ELT(idflist,0,idvec);
	SET_VECTOR_ELT(idflist,1,geneIdVec);
	SET_VECTOR_ELT(idflist,2,transIdVec);

	// Column Names
	SEXP icol_names;
	PROTECT(icol_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(icol_names,0,mkChar("id"));
	SET_STRING_ELT(icol_names,1,mkChar("gene_id"));
	SET_STRING_ELT(icol_names,2,mkChar("transcript_id"));
	setAttrib(idflist,R_NamesSymbol,icol_names);

	// Row Names
	SEXP irow_names;
	PROTECT(irow_names=allocVector(STRSXP,nRows));
	++nProtected;

	for(i=0; i<nRows; ++i)
	{
		sprintf(buf,"%lu",i+1);
		SET_STRING_ELT(irow_names,i,mkChar(buf));
	}
	setAttrib(idflist,R_RowNamesSymbol,irow_names);
	setAttrib(idflist,R_ClassSymbol,mkString("data.frame"));

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Merge id, type, value into adflist data.frame
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

	nCols = 3;
	nRows = (unsigned) l->size;

	// Column 0: id
	SEXP aid_vector;
	PROTECT(aid_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 1: type
	SEXP atype_vector;
	PROTECT(atype_vector=allocVector(STRSXP,nRows));
	++nProtected;

	// Column 2: value
	SEXP aval_vector;
	PROTECT(aval_vector=allocVector(STRSXP,nRows));
	++nProtected;

	ptr_pair_list_rewind(l);
	const ptr_pair_element *e;

	for(i=0;i<nRows;++i)
	{
		e=ptr_pair_list_get_next_element(l);
		INTEGER(aid_vector)[i] = (int) e->id;
		SET_STRING_ELT(atype_vector,i,Rf_mkCharLen(e->first,e->first_len));
		SET_STRING_ELT(aval_vector,i,Rf_mkCharLen(e->second,e->second_len));
	}

	SEXP adflist;
	PROTECT(adflist=allocVector(VECSXP,nCols));
	++nProtected;

	SET_VECTOR_ELT(adflist,0,aid_vector);
	SET_VECTOR_ELT(adflist,1,atype_vector);
	SET_VECTOR_ELT(adflist,2,aval_vector);

	// Column names
	SEXP acol_names;
	PROTECT(acol_names=allocVector(STRSXP,nCols));
	++nProtected;
	SET_STRING_ELT(acol_names,0,mkChar("id"));
	SET_STRING_ELT(acol_names,1,mkChar("type"));
	SET_STRING_ELT(acol_names,2,mkChar("value"));
	setAttrib(adflist,R_NamesSymbol,acol_names);

	// Row names
	SEXP arow_names;
	PROTECT(arow_names=allocVector(STRSXP,nRows));
	++nProtected;
	for(i=0;i<nRows;++i)
	{
		sprintf(buf,"%lu",i+1);
		SET_STRING_ELT(arow_names,i,mkChar(buf));
	}
	setAttrib(adflist,R_RowNamesSymbol,arow_names);

	// Class symbol
	setAttrib(adflist,R_ClassSymbol,mkString("data.frame"));

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// ans result list
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	int list_len=2;
	SEXP ans;
	PROTECT(ans=allocVector(VECSXP,list_len));
	++nProtected;

	SET_VECTOR_ELT(ans,0,idflist);
	SET_VECTOR_ELT(ans,1,adflist);

	// names
	SEXP ans_names;
	PROTECT(ans_names=allocVector(STRSXP,list_len));
	++nProtected;
	SET_STRING_ELT(ans_names,0,mkChar("fixed"));
	SET_STRING_ELT(ans_names,1,mkChar("variable"));
	setAttrib(ans,R_NamesSymbol,ans_names);

	// Class symbol
	setAttrib(ans,R_ClassSymbol,mkString("list"));

	// Return
	ptr_pair_list_destroy(l);
	UNPROTECT(nProtected);
	return ans;
}
*/

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Import data from GTF files
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// First eight columns:
// 1	seqname
// 2	source
// 3	feature	e.g. "CDS", "start_codon", "stop_codon", "exon"
// 4	start 	(1-based)
// 5	end 		(inclusive)
// 6	score	(0 - 1000), No score = "."
// 7	frame	(0-2), Not a coding exon: "."
// 8	group = attributes (variable content)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //


static void print_progress(size_t n) { Rprintf("\r[GTF] %8lu lines processed.", n); }

static void do_fill_char(size_t s, const char* value, void *o)
{
	atmptr<char> *column = (atmptr<char> *) o;
	// Change from 1-based to 0-based index
	column->set( ((int) s) - 1, value);
}

static void do_fill_int(size_t s, const char* value, void *o)
{
	atmptr<int> *column = (atmptr<int> *) o;
	// Change from 1-based to 0-based index
	(*column)[((int) s) - 1] = (int) atol(value);
}


SEXP read_gtf(SEXP pParam, SEXP pProgress)
{
	  if( !Rf_isString(pParam) )
	      error("first argument must be a string, found %s",
	            type2char(TYPEOF(pParam)));

	atmptr<int> pi(pProgress);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Extract pParam values
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	std::ifstream infile(CHAR(STRING_ELT(pParam, 0)));
	if(!infile.is_open())
		error("File not found.");
	char comment = CHAR(STRING_ELT(pParam, 1))[0];
	char delim = CHAR(STRING_ELT(pParam, 2))[0];


	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Open file and process file content
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	gtf::gtf_file gtf(delim);
	gtf.process_lines(infile, pi[0], print_progress, comment);
	Rprintf("\n");
	infile.close();


	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Create empty data.frame for output
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	typedef std::list<std::pair<size_t, std::string> > lps;
	lps map_data = gtf.map_data();
	unsigned nrow = (unsigned) gtf.size();
	unsigned  n_map = (unsigned) map_data.size();
	unsigned ncol = n_map + 9;
	data_frame dfr(nrow, ncol);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Add column data to data.frame
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	dfr.addIdColumn();

	// Add fixed GFF columns
	atmptr<char> seqname(nrow);
	gtf.fill_seqname(do_fill_char, &seqname);
	dfr.addColumn(seqname, "seqid");

	atmptr<char> source(nrow);
	gtf.fill_source(do_fill_char, &source);
	dfr.addColumn(source, "source");

	atmptr<char> feature(nrow);
	gtf.fill_feature(do_fill_char, &feature);
	dfr.addColumn(feature, "feature");

	atmptr<int> start(nrow);
	gtf.fill_start(do_fill_int, &start);
	dfr.addColumn(start, "start");

	atmptr<int> end(nrow);
	gtf.fill_end(do_fill_int, &end);
	dfr.addColumn(end, "end");

	atmptr<char> score(nrow);
	gtf.fill_score(do_fill_char, &score);
	dfr.addColumn(score, "score");

	atmptr<char> strand(nrow);
	gtf.fill_strand(do_fill_char, &strand);
	dfr.addColumn(strand, "strand");

	atmptr<char> frame(nrow);
	gtf.fill_frame(do_fill_char, &frame);
	dfr.addColumn(frame, "frame");

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Add (variable) attribute - columns data to data.frame
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	lps::const_iterator iter;
	for(iter=map_data.begin(); iter != map_data.end(); ++iter)
	{
		atmptr<char> col_data(nrow, true);
		gtf.proc_attr_feature(iter->second.c_str(), do_fill_char, &col_data);
		dfr.addColumn(col_data, iter->second.c_str());
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Return
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	return dfr;
}





SEXP get_cum_max(SEXP pVal)
{
	int i,n;
	int *val, *res;

	// Check type of incoming args
	  if( !Rf_isInteger(pVal) )
	      error("argument must be an integer, found %s",
	            type2char(TYPEOF(pVal)));


	n = LENGTH(pVal);
	if(n == 0)
		return R_NilValue;

	SEXP pRes;
	PROTECT(pRes = allocVector(INTSXP, n));

	val = INTEGER(pVal);
	res = INTEGER(pRes);

	res[0] = val[0];
	for(i = 1; i < n; ++i)
		res[i] = val[i] > res[i-1] ? val[i] : res[i-1];

	UNPROTECT(1);
	return pRes;
}


SEXP gap_overlap(SEXP pQid, SEXP pQlstart, SEXP pQlend, SEXP pQrstart, SEXP pQrend,
		SEXP pRid, SEXP pRlstart, SEXP pRlend, SEXP pRrstart, SEXP pRrend, SEXP pRmaxRend)
{
	/*
	 * Gap-sites:			Represent data on splice-sites:
	 * 						Intronic regions which are flanked by exons.
	 *
	 * Position:			Positions are given as integer values:
	 * 						Coordinates on genomic sequences which
	 * 						ascend from left to right.
	 *
	 * Gap-site records:	Gap-sites are described by four coordinates:
	 * 						lstart	: left  start
	 * 						lend	: left  end
	 * 						rstart	: right start
	 * 						rend	: right end
	 *
	 * 						The coordinates give the genomic positions of
	 * 						the exons which flank the intron.
	 *
	 *
	 * 						lstart        lend     rstart        rend
	 *                           |        |             |        |
	 *                           xxxxxxxxxx             xxxxxxxxxx
	 *                             (exon)     (intron)    (exon)
	 *
	 *
	 *                      Therefore, it is assumed for each record:
	 *                      lstart < lend < rstart < rend
	 *
	 *                      Instead, the RmaxRend value is used, which is
	 *                      the cumulative maximum of all rend values
	 *                      from beginning to the first to the present
	 *                      record.
	 *                      The value can be calculated from (nearly sorded !)
	 *                      rend values by using the previous "get_cum_max"
	 *                      function.
	 *
	 *                      Furthermore, the algorithm only works efficiantly
	 *                      when the records are somehow ordered
	 *                      (e.g. by lstart or lend and rend)
	 *
	 *
	 *                      The algorithm searches for overlaps between
	 *                      query and reference gap-sites as follows:
	 *
	 *						For each query gap-site, the algorithm walks
	 *						reverse from the current position until no more
	 *						overlaps are possible (-> RmaxRend).
	 *
	 *						Then it advances again until the first overlap (hit)
	 *						is identified.
	 *
	 *						From there, the search traverses the reference
	 *                      downstream until no further overlap is possible.
	 *
	 *                      All hits are evaluated by the sum of distances (sod) of
	 *                      the inner boundaries (lend, rstart) between query
	 *                      and reference record.
	 *
	 *                      The best hit is the one with the minimal sod (possibly = 0)
	 *                      which will be returned as result.
	 *
	 *
	 */



	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Check input for expected types
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	if((TYPEOF(pQid) != INTSXP) | (TYPEOF(pQlstart) != INTSXP) | (TYPEOF(pQlend) != INTSXP) |
			(TYPEOF(pQrstart) != INTSXP) | (TYPEOF(pQrend) != INTSXP))
		error("[gap_overlap] Query values must be INTEGER!");

	if( (TYPEOF(pRid) != INTSXP) | (TYPEOF(pRlstart) != INTSXP) | (TYPEOF(pRlend) != INTSXP) |
			(TYPEOF(pRrstart) != INTSXP) | (TYPEOF(pRrend) != INTSXP) | (TYPEOF(pRmaxRend) != INTSXP))
		error("[gap_overlap] Reference values must be INTEGER!");


	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Declaration of variables
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //


	// + + + + + + + + + + + + + + + //
	// 1. Internal values
	// + + + + + + + + + + + + + + + //

	// Running indices
	int i, j;
	// Number of query and reference records
	int nq, nr;
	// Buffer for data.frame names
	char buf[20];

	// Intermediate values for scoring hits in reference records
	int ldiff, rdiff, ldist, rdist, sod;


	// + + + + + + + + + + + + + + + //
	// 2. Input data from arguments
	// + + + + + + + + + + + + + + + //

	// ID and range-gap positions for query and reference
	int *qid, *qLstart, *qLend, *qRstart, *qRend;
	int *rid, *rLstart, *rLend, *rRstart, *rRend, *rMaxRend;


	// + + + + + + + + + + + + + + + //
	// 3. Output data.frame
	// + + + + + + + + + + + + + + + //
	SEXP dflist, col_names, row_names;

	// Column vectors for data.frame (returned as result)
	SEXP p_res_qid, p_res_rid, p_res_ldiff, p_res_rdiff, p_res_sod,
					p_res_nref, p_res_first_refid, p_res_last_refid, p_res_adv, p_res_rev;

	// Integer pointer: point to column arrays of data.frame
	int *res_qid, *res_rid, *res_ldiff, *res_rdiff, *res_nref,
						*res_sod, *res_first_refid, *res_last_refid, *res_adv, *res_rev;



	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Initialize variable values
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	nq = LENGTH(pQid);
	nr = LENGTH(pRid);
	const int nCols=10, nRows=nq;

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// check arguments for equal length
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

	if(LENGTH(pQlstart)!=nq || LENGTH(pQlend)!=nq || LENGTH(pQrstart)!=nq ||
			LENGTH(pQrend)!=nq)
		error("[gap_overlap] Query values must have equal length!");

	if(LENGTH(pRlstart)!=nr || LENGTH(pRlend)!=nr || LENGTH(pRrstart)!=nr ||
			LENGTH(pRmaxRend)!=nr)
		error("[gap_overlap] Reference values must have equal length!");


	// Query value arguments
	qid			= INTEGER(pQid);
	qLstart 	= INTEGER(pQlstart);
	qLend		= INTEGER(pQlend);
	qRstart		= INTEGER(pQrstart);
	qRend		= INTEGER(pQrend);

	// Reference value arguments
	rid			= INTEGER(pRid);
	rLstart		= INTEGER(pRlstart);
	rLend		= INTEGER(pRlend);
	rRstart		= INTEGER(pRrstart);
	rRend		= INTEGER(pRrend);
	rMaxRend	= INTEGER(pRmaxRend);

	// Column vectors for output data.frame
	PROTECT(p_res_qid 			= allocVector(INTSXP, nRows));
	PROTECT(p_res_rid 			= allocVector(INTSXP, nRows));
	PROTECT(p_res_ldiff			= allocVector(INTSXP, nRows));
	PROTECT(p_res_rdiff			= allocVector(INTSXP, nRows));
	PROTECT(p_res_nref 			= allocVector(INTSXP, nRows));
	PROTECT(p_res_sod			= allocVector(INTSXP, nRows));
	PROTECT(p_res_first_refid	= allocVector(INTSXP, nRows));
	PROTECT(p_res_last_refid 	= allocVector(INTSXP, nRows));
	PROTECT(p_res_adv           = allocVector(INTSXP, nRows));
	PROTECT(p_res_rev           = allocVector(INTSXP, nRows));

	// Integer array for output data.frame
	res_qid 		= INTEGER(p_res_qid);
	res_rid 		= INTEGER(p_res_rid);
	res_ldiff		= INTEGER(p_res_ldiff);
	res_rdiff		= INTEGER(p_res_rdiff);
	res_nref 		= INTEGER(p_res_nref);
	res_sod			= INTEGER(p_res_sod);
	res_first_refid	= INTEGER(p_res_first_refid);
	res_last_refid	= INTEGER(p_res_last_refid);
	res_adv			= INTEGER(p_res_adv);
	res_rev			= INTEGER(p_res_rev);



	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Process input data:
	//
	// For each query record: find optimal hit in reference records.
	// i = query index, j = reference index
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	for(i=0, j=0; i<nq; ++i)
	{
		res_qid[i] = qid[i];
		res_rev[i] = 0;
		res_adv[i] = 0;

		// reset to empty values
		//Rprintf("[gap_overlap] New qry: i=%3i qid= %3i j=%3i + + + + + \n", i, res_qid[i], j);

		// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
		// 1. Identify first hit
		// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

		// 1.1 Reference: Step back before first hit
		//Rprintf("[gap_overlap] Reverse: ");
		while( (qLstart[i] <= rMaxRend[j]) && j > 0)
		{
			//Rprintf("%2i ",j);
			--j;
			++res_rev[i];
		}
		//Rprintf("\n");

		// 1.2 Reference: Advance to first hit
		//Rprintf("[gap_overlap] Advance: \n");
		while( j < nr && (qLstart[i] > rMaxRend[j]))
		{
			//Rprintf("[gap_overlap] Advance i: %2i, j: %2i\tqLstart: %5i, rMaxRend: %5i\n",i,j, qLstart[i], rMaxRend[j]);
			++j;
			++res_adv[i];
		}
		//Rprintf("\n");

		// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
		// 2. Past last reference position ?
		//    -> Abort search because no more hits are possible
		// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
		if(j==nr)
		{
			// Write "no-hit" events
			for(;i < nq; ++i)
			{
				res_qid[i]         = qid[i];
				res_rid[i]         = 0;
				res_ldiff[i]       = NA_INTEGER;
				res_rdiff[i]       = NA_INTEGER;
				res_sod[i]         = NA_INTEGER;
				res_nref[i]        = 0;
				res_first_refid[i] = 0;
				res_last_refid[i]  = 0;
				res_adv[i]		   = 0;
				res_rev[i]		   = 0;
			}
			break;
		}

		// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
		// 3. Traverse hit region
		// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

		// 3.1 Set pre-traverse values
		res_rid[i]  =  0; //rid[j];
		res_nref[i] =  0;
		res_sod[i]  =  -1;
		res_first_refid[i] = rid[j];

		// 3.2 Do traverse + identify best hit (= lowest sod)
		//Rprintf("[gap_overlap] i=%2i j=%2i qid=%2i Hits: ", i,j, res_qid[i]);
		while((qRend[i] >= rLstart[j]) && j < nr)
		{
			// Actual overlap ?
			if(qLstart[i] < rRend[j] && qRend[i] > rLstart[j])
			{
				//Rprintf("%2i ",rid[j]);
				++res_nref[i];

				// Differences (may be < 0)
				ldiff = qLend[i] - rLend[j];
				rdiff =  qRstart[i] - rRstart[j];
				// Distance (>= 0)
				ldist = ldiff < 0 ? (-ldiff) : ldiff;
				rdist = rdiff < 0 ? (-rdiff) : rdiff;
				// Sum of distances
				sod   = ldist + rdist;

				// New hit is better than all previous ones
				if( (sod < res_sod[i]) | (res_sod[i] < 0))
				{
					res_rid[i] = rid[j];
					res_ldiff[i] = ldiff;
					res_rdiff[i] = rdiff;
					res_sod[i] = sod;
				}
			}
			++j;
		}
		//Rprintf("\n");


		// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
		// 4. Set post-traverse values
		// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

		// Reverse once for last hit (also: j points not behind array)
		--j;

		// No hits -> Indicate empty result set
		if(res_nref[i] == 0)
		{
			res_sod[i]         = NA_INTEGER; // INT_MIN (currently), see Arith.h
			res_ldiff[i]       = NA_INTEGER;
			res_rdiff[i]       = NA_INTEGER;
			res_first_refid[i] = 0;
			res_last_refid[i]  = 0;
		}
		else
			res_last_refid[i]  = rid[j];
	}


	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Create output data.frame
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

	PROTECT(dflist=allocVector(VECSXP,nCols));

	// Set column vectors
	SET_VECTOR_ELT(dflist, 0, p_res_qid);
	SET_VECTOR_ELT(dflist, 1, p_res_rid);
	SET_VECTOR_ELT(dflist, 2, p_res_ldiff);
	SET_VECTOR_ELT(dflist, 3, p_res_rdiff);
	SET_VECTOR_ELT(dflist, 4, p_res_nref);
	SET_VECTOR_ELT(dflist, 5, p_res_sod);
	SET_VECTOR_ELT(dflist, 6, p_res_first_refid);
	SET_VECTOR_ELT(dflist, 7, p_res_last_refid);
	SET_VECTOR_ELT(dflist, 8, p_res_adv);
	SET_VECTOR_ELT(dflist, 9, p_res_rev);


	// Column Names
	PROTECT(col_names=allocVector(STRSXP,nCols));

	SET_STRING_ELT(col_names, 0, mkChar("qid"));
	SET_STRING_ELT(col_names, 1, mkChar("refid"));
	SET_STRING_ELT(col_names, 2, mkChar("ldiff"));
	SET_STRING_ELT(col_names, 3, mkChar("rdiff"));
	SET_STRING_ELT(col_names, 4, mkChar("nref"));
	SET_STRING_ELT(col_names, 5, mkChar("sod"));
	SET_STRING_ELT(col_names, 6, mkChar("first_refid"));
	SET_STRING_ELT(col_names, 7, mkChar("last_refid"));
	SET_STRING_ELT(col_names, 8, mkChar("nadv"));
	SET_STRING_ELT(col_names, 9, mkChar("nrev"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	// Row Names
	PROTECT(row_names=allocVector(STRSXP,nRows));
	for(i=0;i<nRows;++i)
	{
		sprintf(buf,"%i",i+1);
		SET_STRING_ELT(row_names,i,mkChar(buf));
	}
	setAttrib(dflist,R_RowNamesSymbol,row_names);
	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));

	UNPROTECT(13);
	return dflist;
}


void R_init_refGenome(DllInfo *info)
{
	R_CallMethodDef cmd[] ={
			//{ "split_gtf_attr",			(DL_FUNC) &split_gtf_attr,			2},
			{ "read_gtf",				(DL_FUNC) &read_gtf,                2},
			{ "get_exon_number",		(DL_FUNC) &get_exon_number,       	4},
			{ "get_splice_juncs",		(DL_FUNC) &get_splice_juncs,		4},
			{ "unify_splice_juncs",		(DL_FUNC) &unify_splice_juncs,		9},
			{ "unify_genomic_ranges",	(DL_FUNC) &unify_genomic_ranges,	4},
			{ "overlap_ranges",			(DL_FUNC) &overlap_ranges,			6},
			{ "get_cum_max",			(DL_FUNC) &get_cum_max,				1},
			{ "gap_overlap",			(DL_FUNC) &gap_overlap,			   11},
			{NULL, NULL, 0}
	};
	//			{ "",	(DL_FUNC) &,	}
	R_registerRoutines(info, NULL, cmd, NULL, NULL);
}



} // extern "C"

#endif	/* REFGENOME_CPP_ */
