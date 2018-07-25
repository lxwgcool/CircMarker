#include "clsvelvet.h"
#include "unistd.h"
#include "clsbasealgorithm.h"
#include "clsfastareader.h"
#include <iostream>
#include <stdlib.h>

string ClsVelvet::LocalAssembly(string strFilePath, int iMinOutputLen, int iKmerLen, bool bAppendOrgReads)
{
    if(::access(strFilePath.c_str(), 0) != 0)
    {
        cout << "File is not existed";
        return "";
    }

    //Use Velveth for temperary file generated
    string strCmd = "";
    string strRootPath = GetHigherFolderPath(get_current_dir_name());
    int iKmer = iKmerLen;
    string strVelvetFolderPath = strRootPath + "ThirdPartyTools/Velvet/";
    string strBinVelvetG = strVelvetFolderPath + "velvetg";
    string strBinVelvetH = strVelvetFolderPath + "velveth";
    string strOutput = "Output/LocalAssembly";
    string strLAR = strRootPath + strOutput + "/contigs.fa"; //
    string strFileType = "-fasta";
    string strReadsType = "-short";
    string strScaffoldStatus = "-scaffolding yes"; // make scaffold
    string strMinContigLen = "-min_contig_lgth " + IntToStr(iMinOutputLen);

    //Remove Old contig.fa
    if(::access(strLAR.c_str(), 0) == 0) // 存在
    {
        strCmd = "rm " + strLAR;
        ::system(strCmd.c_str());
    }
    //::system(strCmd.c_str());

    //Remove the old file: do not do it!!! 里面的reads文件是一个很重要的生成后面的velvet的前提，因此不能删掉。
    //strCmd = "rm " + strRootPath + strOutput + "*.*";
    //::system(strCmd.c_str());

    //do velvet
    strCmd = strBinVelvetH + " " + strRootPath + strOutput + " " +
             IntToStr(iKmer) + " " +
             strFileType + " " + strReadsType + " " + strFilePath;
    ::system(strCmd.c_str());
    //Use Velvetg for contig generated
    strCmd = strBinVelvetG + " " +  strRootPath + strOutput + " " +
             strScaffoldStatus + " " + strMinContigLen + " -exp_cov auto";
    ::system(strCmd.c_str());

    //判定一下，看看该文件是否存在，以此防止local assembly不成功的case
    unsigned long ulsize = GetFileSize(strLAR.c_str());
    if(bAppendOrgReads)
    // 这个意味着我们在没生成contigs的时候产生contigs，并且将
    {
        if(ulsize == 0)
        {
            //没生成，使用老板的那个工具去做
            return RelaxMergeByContigsMerger(strFilePath);
        }
        else
        {
            strCmd = "cat " + strFilePath + " >> " + strLAR;
            ::system(strCmd.c_str());
            return strLAR;
        }
    }
    else
    {
        //if(::access(strLAR.c_str(), 0) != 0)
        if(ulsize == 0)
            return "";
        else
            return strLAR;
    }
}

string ClsVelvet::RelaxMergeByContigsMerger(string strFilePath)
{
    unsigned long ulsize = GetFileSize(strFilePath.c_str());
    if(ulsize > 0 && ulsize < 60 * 10 * 1000)
    {
        //ContigsMerger -s 0.2 -i1 -6.0 -i2 -6.0 -x 15 -k 10 -t 15 -m 1 -o merge.info input.fa > output.merged.fa
        //--> Use the Tools ContigsMerger rather than Velvet to Deal with this Problem
        string strCmd = "";
        string strRootPath = GetHigherFolderPath(get_current_dir_name());
        string strOutput = "Output/LocalAssembly/";
        string strLAR = strRootPath + strOutput + "contigs.fa";
        string strContigsMergerPath = strRootPath + "ThirdPartyTools/ContigsMerger";
        //-->Parameter Setting
        float fSilimarRange = 0.2; // -s cutoff for similar regions
        int iOverlapLen = 15;      // -x overlap length
        int iThreadNum = 4;        // -t number of threads
        int iSeedNum = 1;          // -m number of seeds used for quick checking
        int iSeedLen = 10;         // -k seed length
        int iError1 = -6.0;        // -i1 Error1 value
        int iError2 = -6.0;

        // -i2 Error2 value
        strCmd = strContigsMergerPath +
                " -s " + FloatToStr(fSilimarRange, 1) +
                " -i1 " + IntToStr(iError1) +
                " -i2 " + IntToStr(iError2) +
                " -x " + IntToStr(iOverlapLen) +
                " -k " + IntToStr(iSeedLen) +

                " -t " + IntToStr(iThreadNum) +
                " -m " + IntToStr(iSeedNum) +
                " -o merge.info " +
                strFilePath + " > " + strLAR;
        system(strCmd.c_str());

        //-->Delete the Later Part
        /*
        ClsFastaReader* pFastaReader = new ClsFastaReader();
        vector<St_Fasta> vFasta;
        pFastaReader->ReadFastaByKW(strLAR, vFasta, "MERGE");
        delete pFastaReader;
        pFastaReader = NULL;

        strCmd = "rm " + strLAR;
        system(strCmd.c_str());

        ofstream ofs;
        ofs.open(strLAR.c_str());
        for(vector<St_Fasta>::iterator itr = vFasta.begin(); itr != vFasta.end(); itr++)
        {
            ofs << itr->strName << endl;
            ofs << itr->strSeq << endl;
        }
        ofs.close();*/
        return strLAR;
        //<--
    }
    else
        return "";
}
