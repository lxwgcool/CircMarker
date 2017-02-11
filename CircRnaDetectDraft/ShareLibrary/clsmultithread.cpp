#include "clsmultithread.h"
#include "clsbasealgorithm.h"
#include "pthread.h"
#include "ctime"
#include "time.h"
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <map>
#include <stdlib.h>
#include <tr1/memory>
#include <cstdio>
#include <math.h>


/********************St_Super Kmer************************
 * The Functions
 * *******************************************************/
void St_SuperKmer::BuildMapOrg(int iGroupIndex)
{
    mpOrgSimilar.clear();
    if(stInfo.strSeq == "")
        return;
    St_KmerIndex stKmerIndex;
    stKmerIndex.iGroupIndex = iGroupIndex;
    stKmerIndex.iCount = 0;
    for(int iPos = 0; iPos < stInfo.strSeq.length(); iPos++)
    {
        string strPos = stInfo.strSeq.substr(iPos, 1);

        //For "A"
        mpOrgSimilar[stInfo.strSeq.replace(iPos, 1, "A")] = stKmerIndex;
        //For "T"
        mpOrgSimilar[stInfo.strSeq.replace(iPos, 1, "T")] = stKmerIndex;
        //For "G"
        mpOrgSimilar[stInfo.strSeq.replace(iPos, 1, "G")] = stKmerIndex;
        //For "C"
        mpOrgSimilar[stInfo.strSeq.replace(iPos, 1, "C")] = stKmerIndex;

        //Recover the original value
        stInfo.strSeq.replace(iPos, 1, strPos.c_str());
    }
    mpOrgSimilar[stInfo.strSeq].iCount = stInfo.iNum;
}

/*
void St_SuperKmer::BuildMap()
{
    mpSimilar.clear();
    if(stInfo.strSeq == "")
        return;

    for(int iPos = 0; iPos < stInfo.strSeq.length(); iPos++)
    {
        string strPos = stInfo.strSeq.substr(iPos, 1);
        unsigned int iLen = stInfo.strSeq.length();
        //for "A"
        string strCurSimilarSeq = stInfo.strSeq.replace(iPos, 1, "A");
        int iKmerPosStart = 0;
        uint64_t vValue = 0;
        FormKmerTypeShortSeg(strCurSimilarSeq.c_str(), iKmerPosStart, iLen, vValue);
        mpSimilar[vValue] = 0;
        //iKmerPosStart++;

        //For "T"
        strCurSimilarSeq = stInfo.strSeq.replace(iPos, 1, "T");
        FormKmerTypeShortSeg(strCurSimilarSeq.c_str(), iKmerPosStart, iLen, vValue);
        mpSimilar[vValue] = 0;
        //iKmerPosStart++;

        //For "G"
        strCurSimilarSeq = stInfo.strSeq.replace(iPos, 1, "G");
        FormKmerTypeShortSeg(strCurSimilarSeq.c_str(), iKmerPosStart, iLen, vValue);
        mpSimilar[vValue] = 0;
        //iKmerPosStart++;

        //For "C"
        strCurSimilarSeq = stInfo.strSeq.replace(iPos, 1, "C");
        FormKmerTypeShortSeg(strCurSimilarSeq.c_str(), iKmerPosStart, iLen, vValue);
        mpSimilar[vValue] = 0;

        //Back to the original value
        stInfo.strSeq.replace(iPos, 1, strPos.c_str());
    }
}

int St_SuperKmer::GetSize()
{
    int iSize = 0;
    for(map<KmerTypeShort, int>::iterator itr = mpSimilar.begin(); itr != mpSimilar.end(); itr++)
    {
        if(itr->second != 0)
            iSize++;
    }
    return iSize;
}*/

int St_SuperKmer::GetSizeOrg()
{
    int iSize = 0;
    for(map<string, St_KmerIndex>::iterator itr = mpOrgSimilar.begin(); itr != mpOrgSimilar.end(); itr++)
    {
        if(itr->second.iCount != 0)
            iSize += itr->second.iCount;
    }
    if(iSize == 0)
    {
        int i = 0;
        i++;
    }
    return iSize;
}

//we do multi thread here !! -->
std::string exec(const char* cmd)
{
    std::tr1::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while (!feof(pipe.get())) {
        memset(buffer, 0, 128);
        if (fgets(buffer, 128, pipe.get()) != NULL)
        {
            result += buffer;
        }
    }
    return result;
}

void* ClsThreadFunc::SimilarDetect(void* args)
{
    St_SimilarKmer* pStSimilarKmer = (St_SimilarKmer*) args;
    // Search in Regular Container
    int iTemp = 0;
    for(vector<St_KmerInfo>::iterator itrReg = pStSimilarKmer->vRegKmer.begin();
        itrReg != pStSimilarKmer->vRegKmer.end(); itrReg++)
    {
        bool bFind = false;
        if(pStSimilarKmer->mpSuperKmer.find(itrReg->strSeq) != pStSimilarKmer->mpSuperKmer.end()) //find it
        {
            pStSimilarKmer->mpSuperKmer[itrReg->strSeq].iCount += itrReg->iNum;
            bFind = true;
        }
        if(bFind)
        {            
            if(itrReg->strSeq == "ACCTTCATTTA")
            {
                cout << "ACCTTCATTTA" << endl;
            }
            pStSimilarKmer->vRegKmer.erase(itrReg);
            itrReg--;
        }
        iTemp++;
        /*
        for(vector<St_SuperKmer>::iterator itr = mpSuperKmer.begin(); itr != mpSuperKmer.end(); itr++)
        {
            //KmerTypeShort vValue = 0;
            //FormKmerTypeShortSeg(itrReg->strSeq.c_str(), 0, itrReg->strSeq.length(), vValue);
            if(itr->mpOrgSimilar.find(itrReg->strSeq) !=  itr->mpOrgSimilar.end()) //find it
            {
                itr->mpOrgSimilar[itrReg->strSeq]++;
                bFind = true;
                break;
            }
        }
        if(bFind)
        {
            vKmerRegular.erase(itrReg);
            itrReg--;
        }*/
        //iTemp++;
    }
}

ClsMultiThread::ClsMultiThread()
{
}

void ClsMultiThread::FindSimilarKmer(map<string, St_KmerIndex>& mpCombine, vector<St_KmerInfo>& vKmerRegular)
{
    //1: Get the number of CPU
    string strNum = exec("nproc");
    if(strNum.substr(strNum.length()-1,1) == "\n")
        strNum = strNum.substr(0, strNum.length()-1);
    int iCpuNum = atoi(strNum.c_str());

    //2: Divid KmerRegular Table evenly by the number of cpu
    int iLen = ceil((float)vKmerRegular.size() / iCpuNum);
    vector<St_SimilarKmer> vSubGroupKmer;
    vSubGroupKmer.resize(iCpuNum);
    for(int i=0; i<iCpuNum; i++)
    {
        int iStart = i * iLen;
        int iEnd = (i+1)*iLen;
        if(iEnd < vKmerRegular.size())
            vSubGroupKmer[i].vRegKmer.insert(vSubGroupKmer[i].vRegKmer.end(),
                                             vKmerRegular.begin() + iStart, vKmerRegular.begin() + iEnd);
        else
            vSubGroupKmer[i].vRegKmer.insert(vSubGroupKmer[i].vRegKmer.end(),
                                             vKmerRegular.begin() + iStart,
                                             vKmerRegular.end());
        vSubGroupKmer[i].mpSuperKmer.insert(mpCombine.begin(), mpCombine.end());
    }

    //Create Threads
    int start_s = time(NULL);
    pthread_t tids[iCpuNum];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for(int i=0; i<iCpuNum; i++)
    {
        int ret = pthread_create(&tids[i], &attr, ClsThreadFunc::SimilarDetect, (void*)&(vSubGroupKmer[i]));
        if(ret != 0)
            cout << "Thread error: " << ret << endl;
    }
    pthread_attr_destroy(&attr);

    void* status;
    for(int i=0; i<iCpuNum; i++)
    {
        int ret = pthread_join(tids[i], &status);
        if(ret != 0)
            cout << "Thread Error: " << ret << endl;
        else
            cout << "Thread Status: " << (long)status << endl;
    }
    //pthread_exit(NULL);
    int stop_s = time(NULL);
    double dRuningTime = difftime(stop_s, start_s);
    string strTimeFormat = ::GetHMSTimeFormat(dRuningTime);
    cout << "Running Time: " << strTimeFormat << endl;

    //Arrange the kmer
    for(map<string, St_KmerIndex>::iterator itr = mpCombine.begin(); itr != mpCombine.end(); itr++)
    {
        int iBase = itr->second.iCount;
        int iAccumulate = 0;
        for(vector<St_SimilarKmer>::iterator subItr = vSubGroupKmer.begin();
            subItr != vSubGroupKmer.end(); subItr++)
        {
            iAccumulate += subItr->mpSuperKmer[itr->first].iCount - iBase;
        }
        itr->second.iCount = iBase + iAccumulate;
    }
    //好折腾啊！！！！！！！！
    //2: 将剩下的没能归组的再进行合并，得到总的没有归组的kmer
    vKmerRegular.clear();
    for(vector<St_SimilarKmer>::iterator itr = vSubGroupKmer.begin();
        itr != vSubGroupKmer.end(); itr++)
    {
        vKmerRegular.insert(vKmerRegular.end(), itr->vRegKmer.begin(), itr->vRegKmer.end());
    }
}
