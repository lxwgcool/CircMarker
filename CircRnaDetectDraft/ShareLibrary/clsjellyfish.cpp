#include "clsjellyfish.h"
#include "clsbasealgorithm.h"
#include "unistd.h"
#include "clsfastareader.h"
#include <iostream>
#include <stdlib.h>

string ClsJellyfish::GetSolidKmer(string strReadsFaPath, int iKmerLen, bool bForceNew,
                                  string strAnchorKmerFileName,
                                  bool bUseFixThreshould, int iCoverageThreshold)
{
    string strRoot = "../../ShareLibrary/"; //GetHigherFolderPath(GetCurExeFolderPath());

    //Target is get the AnchorKmer.fa
    //strRoot + "ThirdPartyTools/Jellyfish/data/" + strAnchorKmerFileName;
    string strAnchorKmerFasta = strRoot + "Jellyfish/data/" + strAnchorKmerFileName;
    if(!bForceNew && access(strAnchorKmerFasta.c_str(), 0) == 0)
    {
        cout << "Anchor Kmer Exsited!" << endl;
        return strAnchorKmerFasta;
    }

    //1: Run JellyFish, create mer_count.jf file -->
    string strJellyFishPath = strRoot + "Jellyfish/bin/jellyfish";
    string strKmerLen = IntToStr(iKmerLen);
    string strThreadsNum = IntToStr(4);
    string strHushElementNum = "150M";
    string strMerCountFile = strRoot + "Jellyfish/data/mer_counts.jf";


    string strCmd = strJellyFishPath + " count -m " + strKmerLen + " -s " + strHushElementNum +
                    " -o " + strMerCountFile +
                    " -t " + strThreadsNum + " -C " + strReadsFaPath;
    if(bForceNew || access(strMerCountFile.c_str(), 0) != 0) // do not exist
        system(strCmd.c_str());

    //2: dump the kmer information out of mer_count.jf
    string strDumpFile = strRoot + "Jellyfish/data/mer_counts_dumps.fa";
    strCmd = strJellyFishPath + " dump " + strMerCountFile +
             " > " + strDumpFile;
    if(bForceNew || access(strDumpFile.c_str(), 0) != 0) // do not exist
        system(strCmd.c_str());

    //3: Create Histogram
    string strHistoPath = strRoot + "Jellyfish/data/histo.ini";
    strCmd = strJellyFishPath + " histo " + strMerCountFile + " > " + strHistoPath;
    system(strCmd.c_str());

    //Read Histogram file
    ifstream ifs;
    ifs.open(strHistoPath.c_str());
    string strLine;
    int iSumFreq = 0;
    int iSumKmerNum = 0;
    getline(ifs, strLine);
    while(!ifs.eof())
    {
        getline(ifs, strLine);
        string::size_type stpzPos = strLine.find(" ");
        string strFreq = strLine.substr(0, stpzPos);
        string strKmerLen = strLine.substr(stpzPos + 1, strLine.length() - stpzPos);
        iSumKmerNum += atoi(strKmerLen.c_str());
        iSumFreq += atoi(strFreq.c_str()) * atoi(strKmerLen.c_str());
    }
    float fAverage = (float)iSumFreq / iSumKmerNum;
    cout << "Average Coverage: " << FloatToStr(fAverage);
    ifs.close();

    //4: Load Fa File and Calculate the sum of kmer and the sum of appeared times.
    ClsFastaReader* pFastaReader = new ClsFastaReader();
    vector<St_Fasta> vKmerFasta;
    pFastaReader->ReadFastaRegular(strDumpFile, vKmerFasta); // collect the kmer which frequency larger than fAverage
    cout << IntToStr(vKmerFasta.size()) << endl;
    vector<St_Fasta> vHighCovKmerFasta;
    for(vector<St_Fasta>::iterator itr = vKmerFasta.begin(); itr != vKmerFasta.end(); itr++)
    {
        if(itr->strSeq == "")
            continue;
        if(atoi(itr->strName.c_str()) > fAverage)
            vHighCovKmerFasta.push_back(*itr);
    }
    vKmerFasta.clear();
    cout << IntToStr(vHighCovKmerFasta.size()) << endl;
    delete pFastaReader;
    pFastaReader = NULL;

    //Save those new Anchor Kmer as a new fasta file
    ofstream ofsAnchorKmer;
    ofsAnchorKmer.open(strAnchorKmerFasta.c_str());
    for(vector<St_Fasta>::iterator itr = vHighCovKmerFasta.begin();
        itr != vHighCovKmerFasta.end(); itr++)
    {
        ofsAnchorKmer << ">" << itr->strName << endl;
        ofsAnchorKmer << itr->strSeq << endl;
    }
    ofsAnchorKmer.close();
    vHighCovKmerFasta.clear();
    return strAnchorKmerFasta;
}

string ClsJellyfish::GetAllKmer(string strReadsFaPath, int iKmerLen)
{
    string strRoot = "../../ShareLibrary/"; //GetHigherFolderPath(GetCurExeFolderPath());

    //1: Run JellyFish, create mer_count.jf file -->
    string strJellyFishPath = strRoot + "Jellyfish/bin/jellyfish";
    string strKmerLen = IntToStr(iKmerLen);
    string strThreadsNum = IntToStr(4);
    string strHushElementNum = "150M";
    string strMerCountFile = strRoot + "Jellyfish/data/mer_counts.jf";


    string strCmd = strJellyFishPath + " count -m " + strKmerLen + " -s " + strHushElementNum +
                    " -o " + strMerCountFile +
                    " -t " + strThreadsNum + " -C " + strReadsFaPath;
    system(strCmd.c_str());

    //2: dump the kmer information out of mer_count.jf
    string strDumpFile = strRoot + "Jellyfish/data/mer_counts_dumps.fa";
    strCmd = strJellyFishPath + " dump " + strMerCountFile +
             " > " + strDumpFile;
    system(strCmd.c_str());

    return strDumpFile;
}
