#ifndef CLSFASTQREADER_H
#define CLSFASTQREADER_H

#include <string>
#include <vector>
using namespace std;

struct St_Fastq
{
    string strName;
    string strComments;
    string strSeq;
    string strQuality;

    St_Fastq()
    {
        Init();
    }

    void Init()
    {
        strName = "";
        strComments = "+";
        strSeq = "";
        strQuality = "";
    }
};

class ClsFastqReader
{
public:
    ClsFastqReader();

public:
    void ReadFastqFile(string strPath, vector<St_Fastq>& vFastq, bool bToUpperCase=true);
    void FastqToFasta(string strFastqFilePath, string strFastaFilePath);
};

#endif // CLSFASTQREADER_H
