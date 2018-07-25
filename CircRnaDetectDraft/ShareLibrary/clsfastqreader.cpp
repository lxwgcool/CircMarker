#include "clsfastqreader.h"
#include <zlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "../ShareLibrary/clsbasealgorithm.h"
#include "../ShareLibrary/FastqFileParse/kseq.h"
KSEQ_INIT(gzFile, gzread)

ClsFastqReader::ClsFastqReader()
{
}

void ClsFastqReader::ReadFastqFile(string strPath, vector<St_Fastq>& vFastq, bool bToUpperCase)
{
    //vFastq.clear();
    St_Fastq stFastq;

    gzFile fp;
    kseq_t *seq;
    int l;

    if(strPath == "")
        return;

    if(access(strPath.c_str(), 0) != 0) //File do not existed
    {
        cout << "Error: File does not existed!" << endl;
        return;
    }

    fp = gzopen(strPath.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0)
    {
        stFastq.Init();
        //Record Name
        stFastq.strName = seq->name.s;
        //Record Comments
        if (seq->comment.l)
            stFastq.strComments = seq->comment.s;
        //Record Sequence        
        stFastq.strSeq = seq->seq.s;
        if(bToUpperCase)
        {
            ToUpper(stFastq.strSeq);
        }
        //Record Quality
        if (seq->qual.l)
            stFastq.strQuality = seq->qual.s;
        vFastq.push_back(stFastq);
    }
    //cout << IntToStr(vFastq.size()) << endl;
    kseq_destroy(seq);
    gzclose(fp);
    return;
}

void ClsFastqReader::FastqToFasta(string strFastqFilePath, string strFastaFilePath)
{
    vector<St_Fastq> vFastq;
    //1: Read Fastq File
    ReadFastqFile(strFastqFilePath, vFastq);

    //2: Save Name and sequence value as Fasta file
    ofstream ofs;
    ofs.open(strFastaFilePath.c_str());
    for(vector<St_Fastq>::iterator itr = vFastq.begin(); itr != vFastq.end(); itr++)
    {
        ofs << ">" << itr->strName << endl;
        ofs << itr->strSeq << endl;
    }
    vFastq.clear();
    ofs.close();
}
