#ifndef CLSFINDCANDIDATE_H
#define CLSFINDCANDIDATE_H
#include "../../ShareLibrary/clsgtfparse.h"
#include "../../ShareLibrary/clsfastqreader.h"
#include "clskmertable.h"

enum En_CircType{ctSelf, ctRegular, ctMax};

struct St_Candidate
{
    unsigned char ucChromIndex;
    int iStartPos;
    int iEndPos;
    //-->
    int iSupportNum; // Especially For My result
    //<--
    bool bRC;
    En_CircType enCircType;
    string strTag;

    St_Candidate()
    {
        ucChromIndex = 0;
        iStartPos = 0;
        iEndPos = 0;
        iSupportNum = 0;
        bRC = false;
        enCircType = ctMax;
        strTag = "xx";
    }

    St_Candidate(unsigned char V1, int V2, int V3, int V4=1, bool V5=false,
                 En_CircType V6=ctMax, string V7="")
    {
        ucChromIndex = V1;
        iStartPos = V2;
        iEndPos = V3;
        iSupportNum = V4;
        bRC = V5;
        enCircType = V6;
        strTag = "xx";
    }

    void SetChromStartEnd(unsigned char V1, int V2, int V3, int V4=1)
    {
        this->ucChromIndex = V1;
        this->iStartPos = V2;
        this->iEndPos = V3;
        this->iSupportNum = V4;
    }

    void SetCircType(En_CircType enType)
    {
        this->enCircType = enType;
    }

    string GetTypeString()
    {
        if(enCircType == ctMax)
            return "M"; //Max or Mixture
        else if(enCircType == ctSelf)
            return "S"; //Single
        else if(enCircType == ctRegular)
            return "R"; //Regular
        else
            return "N"; //Nil
    }

    bool operator == (const St_Candidate& rhs) const // 我们在这里不对direction进行比较
    {
        if( this->ucChromIndex == rhs.ucChromIndex &&
            this->iStartPos == rhs.iStartPos &&
            this->iEndPos == rhs.iEndPos )
            return true;
        else
            return false;
    }

    bool operator < (const St_Candidate& rhs) const
    {
        return this->iStartPos < rhs.iStartPos;
        /*
        if(this->ucChromIndex < rhs.ucChromIndex)
            return true;
        else // 如果大于或者等于
        {
            if(this->iStartPos < rhs.iStartPos)
                return true;
            else
            {
                if(this->iEndPos < rhs.iEndPos)
                    return true;
                else
                    return false;
            }
        }*/
    }
};

//---------->
struct St_HitExon
{
    St_PosInfo stPI;
    int iCount;
    int iLen; // we need to know the length
    vector<char> vPart;

    St_HitExon():iCount(0), iLen(-1)
    {}

    St_HitExon(St_PosInfo& stV1, int iV2, int iV3)
    {
        stPI = stV1;
        iCount = iV2;
        iLen = iV3;
        vPart.clear();
        vPart.push_back(stV1.cPart);
    }

    void UpdateVPart(char cPart)
    {
        vPart.push_back(cPart);
    }
};

struct St_HitCase
{
    vector<St_HitExon> vHitExons;

    St_HitCase(){}

    St_HitCase(St_PosInfo& stV1, int iV2, int iV3)
    {
        vHitExons.push_back(St_HitExon(stV1, iV2, iV3));
    }

    int GetHitSum() //这样是对的，因为你不知道是hit了一个，两个还是三个，这都有可能
    {
        int iHitSum = 0;
        for(vector<St_HitExon>::iterator itr = vHitExons.begin(); itr != vHitExons.end(); itr++)
        {
            iHitSum += itr->iCount;
        }
        return iHitSum;
    }

    int GetHitLen() //这样是对的，因为你不知道是hit了一个，两个还是三个，这都有可能
    {
        int iHitLen = 0;
        for(vector<St_HitExon>::iterator itr = vHitExons.begin(); itr != vHitExons.end(); itr++)
        {
            iHitLen += itr->iLen;
        }
        return iHitLen;
    }
};
//<---------

class ClsFindCandidate
{
public:
    ClsFindCandidate();

public:
    void CheckHitting(string strReads1Path, string strReads2Path,
                      int iMinSupportReads, float fKmerRatio,
                      map<unsigned int, vector<St_PosInfo> >& mpKT,
                      vector<St_Row_Chrom>& vChrom); // Check Hitting

private:
    bool CheckSampling(unsigned int* arrySamplingKmer,
                       map<unsigned int, vector<St_PosInfo> >& mpKT,
                       vector<St_Row_Chrom>& vChrom);

    void CheckHitForCurReads(string& strSeq, bool bSeqRC,
                             map<unsigned int, vector<St_PosInfo> >& mpKT,
                             vector<St_Row_Chrom>& vChrom, float fKmerRatio);
    bool CheckSelfCircRNA(St_HitCase* pHitCase, int iReadLen, vector<St_Row_Chrom>& vChrom);
    bool CheckRegularCircRNA(St_HitCase* pHitCase, int iReadLen, vector<St_Row_Chrom>& vChrom);

    void RefineCandiate(int iMinSupportReads);

    void AssembleReads(string strReads1Path, string strReads2Path, vector<St_Fastq>& vFastq);

    bool CheckHitNum(string& strSeq, float fKmerRatio, St_HitCase* pHitCase);
    bool CheckHitPart(vector<St_Row_Chrom>& vChrom, St_HitCase* pHitCase,
                      string& strSeq, float fKmerRatio, bool bSeqRC);

private:
    void UpdateCurCandiForRegularDirection(St_HitCase* pHitCase,
                                           vector<St_Row_Chrom>& vChrom,
                                           St_Candidate& stCurCandi);

private:
    //map<St_Candidate, int> m_mpSelfCircCandi;    // int is the number of support reads
    //map<St_Candidate, int> m_mpRegCircCandi; // int is the number of support reads

    //Use vector to check if map works fine
    vector<St_Candidate> m_vSelfCircCandi;
    vector<St_Candidate> m_vRegCircCandi;

private:
    int m_iTotalNum;
    int m_iRegNum;

};

#endif // CLSFINDCANDIDATE_H
