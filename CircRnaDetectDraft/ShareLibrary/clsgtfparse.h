#ifndef CLSGTFPARSE_H
#define CLSGTFPARSE_H
#include <vector>
#include <string>
using namespace std;

/*
struct St_CandiAtom //Circular RNA Candidate Atom
{
    int iStart;
    int iEnd;

    St_CandiAtom():iStart(-1), iEnd(-1)
    {}

    St_CandiAtom(int iV1, int iV2)
    {
        iStart = iV1;
        iEnd = iV2;
    }

    void Refresh()
    {
        iStart = -1;
        iEnd = -1;
    }
};*/

struct St_Raw_Exon
{
    int iStart;
    int iEnd;
    bool bRC;

    string strHead2Bp;
    string strTail2Bp;

    St_Raw_Exon()
    {
        iStart = -1;
        iEnd = -1;
        bRC = false;
        strHead2Bp = "";
        strTail2Bp = "";
    }

    void Refresh()
    {
        iStart = -1;
        iEnd = -1;
        bRC = false;
        strHead2Bp = "";
        strTail2Bp = "";
    }

    bool GetIsSupportCircRNA()
    {
        bool bSupport = false;
        if(bRC) //反向
        {
            if(strHead2Bp == "AC" || strTail2Bp == "CT")
            {
                bSupport = true;
            }
        }
        else // 正向
        {
            if(strHead2Bp == "AG" || strTail2Bp == "GT")
            {
                bSupport = true;
            }
        }
        return bSupport;
    }

    bool GetIsCRNAHeadExon()
    {
        bool bSupport = false;
        if(bRC) // 反向
        {
            if(strHead2Bp == "AC")
                bSupport = true;
        }
        else //正向
        {
            if(strHead2Bp == "AG")
                bSupport = true;
        }
        return bSupport;
    }

    bool GetIsCRNATailExon()
    {
        bool bSupport = false;
        if(bRC) // 反向
        {
            if(strTail2Bp == "CT")
                bSupport = true;
        }
        else //正向
        {
            if(strTail2Bp == "GT")
                bSupport = true;
        }
        return bSupport;
    }

    int GetLength()
    {
        return iEnd - iStart + 1;
    }
};

struct St_Raw_Transcript
{
    int iStart;
    int iEnd;
    string strID;
    bool bRC;
    vector<St_Raw_Exon> vRExon; //Raw Exon
    //vector<St_CandiAtom> vPossibleCRNA; // possible circular rna

    St_Raw_Transcript()
    {
        iStart = -1;
        iEnd = -1;
        strID = "";
        bRC = false;
    }

    void Refresh()
    {
        iStart = -1;
        iEnd = -1;
        strID = "";
        bRC = false;
        vRExon.clear();
        //vPossibleCRNA.clear();
    }
};

struct St_Raw_Gene
{
    string strID;
    string strName;
    string strChromoson;
    int iStart;
    int iEnd;
    bool bRC;
    vector<St_Raw_Transcript> vRT; // Raw Transcript

    St_Raw_Gene()
    {
        strName = "";
        strChromoson = "";
        iStart = -1;
        iEnd = -1;
        bRC = false;
    }

    void Refresh()
    {
        strName = "";
        strChromoson = "";
        iStart = -1;
        iEnd = -1;
        bRC = false;
        vRT.clear();
    }
};

struct St_Row_Chrom
{
    string strName;
    vector<St_Raw_Gene> vRG; //Row Gene

    St_Row_Chrom(): strName("")
    {}

    void Refresh()
    {
        strName = "";
        vRG.clear();
    }
};

class ClsGTFParse
{
public:
    ClsGTFParse();

public:
    bool ReadGTF(string strFilePath, vector<St_Row_Chrom>& vChrom);

    void GetRNARef(string strDNARef, string strRNARefPath, vector<St_Raw_Gene>& vGenes, bool bExportSeq=true);

    void LoadRNARef(string strRNARef, vector<St_Raw_Gene>& vGenes, bool bLoadSeq=false);

    void GetTagValue(string strDNARef, int iReadsLen, float fKmerRatio, int iKmerLen,
                     vector<St_Row_Chrom>& vChrom);

    void ColletcPossibleCRNA(vector<St_Raw_Gene>& vGenes);
};

#endif // CLSGTFPARSE_H
