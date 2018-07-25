#ifndef CLSKMERTABLE_H
#define CLSKMERTABLE_H

#include "stdint.h"
#include <map>
#include "../../ShareLibrary/clsgtfparse.h"
#include "../../ShareLibrary/clsfastareader.h"

const int KMERLEN = 15; // do not allowed to extend than 16 (we use int_32 to save this value)

struct St_PosInfo  // Notice: 在这里我们不需要记录相应的位置，只需要记录是哪个index就可以了  --> 并且我们不需要记录相应的rc情况，因为我们通过查询可以立马得到
{
    unsigned char ucChromIndex;      // For Human: maximum 256
    uint16_t uiGeneIndex;              // For Human: maximum less than 32768
    unsigned char ucTranscriptIndex; // For Human: maximum should less than 256
    unsigned char ucExonIndex;       // For Human: maximum less than 256
    char cPart;                      // U: unknown, S: Start, E: End 知道该位置是头还是尾

    St_PosInfo()
    {
        ucChromIndex = 0;
        uiGeneIndex = 0;
        ucTranscriptIndex = 0;
        ucExonIndex = 0;
        cPart = 'U'; // Unknown
    }

    St_PosInfo(unsigned char ucV1, uint16_t uiV2, unsigned char ucV3, unsigned char ucV4)
    {
        ucChromIndex = ucV1;
        uiGeneIndex = uiV2;
        ucTranscriptIndex = ucV3;
        ucExonIndex = ucV4;
    }

    bool CheckSameTranscript(const St_PosInfo& rhs)
    {
        if( this->ucChromIndex == rhs.ucChromIndex &&
            this->uiGeneIndex == rhs.uiGeneIndex &&
            this->ucTranscriptIndex == rhs.ucTranscriptIndex )
            return true;
        else
            return false;
    }

    bool CheckSamePart(const St_PosInfo& rhs)
    {
        if(rhs.cPart == this->cPart)
            return true;
        else
            return false;
    }

    bool operator == (const St_PosInfo& rhs) // do not consider cPart
    {
        if( this->ucChromIndex == rhs.ucChromIndex &&
            this->uiGeneIndex == rhs.uiGeneIndex &&
            this->ucTranscriptIndex == rhs.ucTranscriptIndex &&
            this->ucExonIndex == rhs.ucExonIndex )
            return true;
        else
            return false;
    }

    bool operator < (const St_PosInfo& rhs) const    
    {      
        //1: Compare chromosone Index
        if(this->ucChromIndex < rhs.ucChromIndex)
            return true;
        else if (this->ucChromIndex > rhs.ucChromIndex)
            return false;
        else // equal
        {
            //2: Compare Gene Index
            if(this->uiGeneIndex < rhs.uiGeneIndex)
                return true;
            else if (this->uiGeneIndex > rhs.uiGeneIndex)
                return false;
            else // equal
            {
                //3: Compare Transcript Index
                if(this->ucTranscriptIndex < rhs.ucTranscriptIndex)
                    return true;
                else if (this->ucTranscriptIndex > rhs.ucTranscriptIndex)
                    return false;
                else // equal
                {
                    //4: Compare Exon Index
                    if(this->ucExonIndex < rhs.ucExonIndex)
                        return true;
                    else if (this->ucExonIndex > rhs.ucExonIndex)
                        return false;
                    else // equal
                    {
                        return false; // stable sort
                    }
                }
            }
        }
    }
};

class ClsKmerTable
{
public:
    ClsKmerTable();
    ~ClsKmerTable();

public:
    void CreateKmerTable(map<unsigned int, vector<St_PosInfo> >& mpKT,
                         string strDNARef, int iReadsLen, float fKmerRatio,
                         St_Row_Chrom* pRowChrom, St_Fasta* pRef, int iChromIndex);

//    map<unsigned int, vector<St_PosInfo> >& GetKT();

private:
    void CollectKmerInfo( map<unsigned int, vector<St_PosInfo> >& mpKT,
                          St_Fasta* pCurRefFa, St_Raw_Exon* pExon, St_PosInfo& stKmerInfo, int iBoundLen);

private:
    void InsertCurKmer(map<unsigned int, vector<St_PosInfo> >& mpKT,
                       string strCurKmer, St_PosInfo& stPosInfo);
//    map<unsigned int, vector<St_PosInfo> > m_mpKT;
};

#endif // CLSKMERTABLE_H
