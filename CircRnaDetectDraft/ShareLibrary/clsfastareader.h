#ifndef CLSFASTAREADER_H
#define CLSFASTAREADER_H
#include <vector>
#include <string>
using namespace std;

enum En_SegmentType{stNone=0, stContig, stScaf, stMax};
enum En_FastaType{ftNone=0, ftContig, ftScaf, ftQuery, ftRef, ftMax};
struct St_Fasta
{
    string strName;
    string strSeq;
    En_SegmentType enType;

    St_Fasta():strName(""),strSeq(""), enType(stNone)
    {
    }

    St_Fasta(string strNameValue, string strSeqValue)
    {
        strName = strNameValue;
        strSeq = strSeqValue;
        enType = stNone;
    }

    void Init()
    {
        strName = "";
        strSeq = "";
        enType = stNone;
    }

    bool operator == (St_Fasta& stCmp)
    {
        if(this->strName == stCmp.strName &&
           this->strSeq == stCmp.strSeq )
            return true;
        else
            return false;
    }

    ~St_Fasta(){}
};

class ClsFastaReader
{
public:
    ClsFastaReader();
    ~ClsFastaReader();
public:
    int ReadFastaRegular(string& strPath, vector<St_Fasta>& vFasta, bool bRecordSeq=true); // normal read fasta
    // read fast file by keywards
    int ReadFastaByKW(string& strPath, vector<St_Fasta>& vFasta, string strKeyWord);
    // read fast total
    int ReadFastaDraftGeno(string& strPath,
                           vector<St_Fasta>& vFasta,
                           bool bRecordSeq = true); // Contain both scaffold and contigs
};

#endif // CLSFASTAREADER_H
