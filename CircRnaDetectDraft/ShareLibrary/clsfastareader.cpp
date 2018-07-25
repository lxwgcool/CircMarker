#include "clsfastareader.h"
#include "clsbasealgorithm.h"
#include "unistd.h" // For API: "access"
#include <iostream>
#include <fstream>
#include <string.h>

const char* CZSCAF = "scaffold"; //scaffold
const char* CZCNTG = "contig"; //contig


ClsFastaReader::ClsFastaReader()
{
}

ClsFastaReader::~ClsFastaReader()
{
}

// normal read fasta --> 读取所有的存储在fasta文件中的序列信息
// bRecordSeq 意味着要不要把fa中的序列记录到结构体中(因为有的时候，我们只需要名字不需要记录序列信息，而且序列信息size较大)
int ClsFastaReader::ReadFastaRegular(string& strPath, vector<St_Fasta>& vFasta, bool bRecordSeq)
{
    if(::access(strPath.c_str(), 0) != 0)
    {
        cout << "Error:Fasta File do not existed!" << endl;
        return 1;
    }
    vFasta.clear();
    //Parse the contig file
    ifstream infile;
    infile.open(strPath.c_str(), ios::in);
    string strLine = "";
    St_Fasta stFasta;
    while(!infile.eof()) // check if reached the end
    {
        getline(infile, strLine);
        //infile >> strLine;
        //if(strLine.find('.') != string::npos)
        //    continue;
        string::size_type sztp = strLine.find('>', 0);
        if(sztp != string::npos) // this line is the name of contig
        {
            if(stFasta.strName != "")
            {
                vFasta.push_back(stFasta);
                stFasta.strName = "";
                stFasta.strSeq = "";
            }
            //Set Name
            stFasta.strName = strLine.substr(1, strLine.length()-1);
            stFasta.strSeq = "";
            continue;
        }
        //Set Sequence
        if(infile.eof()) // for the last item
        {
            if(stFasta.strName != "")
                vFasta.push_back(stFasta);
        }
        else
        {
            if(bRecordSeq)
            {
                //convert every string to Upper Case
                ToUpper(strLine);
                stFasta.strSeq += strLine;
            }
        }
    }
    infile.close();
    return 0;
}

// read fast file by keywards，仅仅只记录跟keywords能够项符合的条目
int ClsFastaReader::ReadFastaByKW(string& strPath, vector<St_Fasta>& vFasta, string strKeyWord)
{
    if(::access(strPath.c_str(), 0) != 0)
    {
        cout << "Error:Fasta File do not existed!" << endl;
        return 1;
    }
    vFasta.clear();
    //Parse the contig file
    ifstream infile;
    infile.open(strPath.c_str(), ios::in);
    string strLine = "";
    St_Fasta stFasta;
    while(!infile.eof()) // check if reached the end
    {
        getline(infile, strLine);
        //infile >> strLine;
        //if(strLine.find('.') != string::npos)
        //    continue;
        string::size_type sztp = strLine.find('>', 0);
        if(sztp != string::npos) // this line is the name of contig
        {
            if(stFasta.strName != "" &&
               stFasta.strName.find(strKeyWord.c_str(), 0) != string::npos ) // 只有后面有scaffold字样的才是
            {
                if(strcmp(strKeyWord.c_str(), CZSCAF) == 0)
                    stFasta.enType = stScaf;
                else if(strcmp(strKeyWord.c_str(), CZCNTG) == 0)// 默认不是scaffold的片段，那么就都是contig
                    stFasta.enType = stContig;
                else
                    stFasta.enType = stNone;
                vFasta.push_back(stFasta);

                //初始化fasta file的值
                stFasta.Init();
            }
            //Set Name
            stFasta.strName = strLine.substr(1, strLine.length()-1);
            stFasta.strSeq = "";
            continue;
        }
        //Set Sequence
        if(infile.eof()) // for the last item
        {
            if(stFasta.strName.find(strKeyWord.c_str(), 0) != string::npos)
            {
                if(strcmp(strKeyWord.c_str(), CZSCAF) == 0)
                    stFasta.enType = stScaf;
                else if(strcmp(strKeyWord.c_str(), CZCNTG) == 0)// 默认不是scaffold的片段，那么就都是contig
                    stFasta.enType = stContig;
                else
                    stFasta.enType = stNone;
                vFasta.push_back(stFasta);
            }
        }
        else
        {
            stFasta.strSeq += strLine;
        }
    }
    infile.close();
    return 0;
}

// read fast total
// 对于 DraftGeno，不是scaffold就是contig
int ClsFastaReader::ReadFastaDraftGeno(string& strPath,
                                       vector<St_Fasta>& vFasta,
                                       bool bRecordSeq) // Contain both scaffold and contigs
{
    if(::access(strPath.c_str(), 0) != 0)
    {
        cout << "Error:Fasta File do not existed!" << endl;
        return 1;
    }
    vFasta.clear();
    //Parse the contig file
    ifstream infile;
    infile.open(strPath.c_str(), ios::in);
    string strLine = "";
    St_Fasta stFasta;
    while(!infile.eof()) // check if reached the end
    {
        getline(infile, strLine);
        //if(strLine.find('.') != string::npos)
        //    continue;
        string::size_type sztp = strLine.find('>', 0);
        if(sztp != string::npos) // this line is the name of contig
        {
            if(stFasta.strName != "") // 只有后面有scaffold字样的才是
            {
                //-->Record current Fasta Unit
                if(stFasta.strName.find(CZSCAF) != string::npos)
                    stFasta.enType = stScaf;
                else // 默认不是scaffold的片段，那么就都是contig
                    stFasta.enType = stContig;
                vFasta.push_back(stFasta);
                //<--

                stFasta.strName = "";
                stFasta.strSeq = "";
            }
            //Set Name
            stFasta.strName = strLine.substr(1, strLine.length()-1);
            if(stFasta.strName.find(" ") != string::npos) // 存在空格
            {
                stFasta.strName = stFasta.strName.substr(0, stFasta.strName.find(" "));
            }
            stFasta.strSeq = "";
            continue;
        }
        //Set Sequence
        if(infile.eof()) // for the last item
        {
            if(stFasta.strName != "")
            {
                //-->Record current Fasta Unit
                if(stFasta.strName.find(CZSCAF, 0) != string::npos)
                    stFasta.enType = stScaf;
                else // 默认不是scaffold的片段，那么就都是contig
                    stFasta.enType = stContig;
                vFasta.push_back(stFasta);
            }
        }
        else
        {
            if(bRecordSeq)
            {
                stFasta.strSeq += strLine;
            }
        }
    }
    infile.close();
    return 0;
}

