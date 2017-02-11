/*
 * Author: lxwg
 * Notice: you need to add the line below in your project configuration file:
 * ----------->
 * CONFIG += thread #This is to support Multi_Thread (pthread!!!)
 * <-----------
 */

#ifndef CLSMULTITHREAD_H
#define CLSMULTITHREAD_H
#include <map>
#include <string>
//#include "KmerUtils.h"
#include <vector>
#include <map>
#include <stdint.h>
using namespace std;

struct St_KmerInfo
{
    string strSeq;
    unsigned int iID;
    int iNum; // 这个就是频率 ！！！
    St_KmerInfo():strSeq(""), iID(0), iNum(0)
    {}
};

struct St_KmerIndex
{
    int iGroupIndex;
    int iCount;

    St_KmerIndex():iGroupIndex(-1), iCount(-1)
    {}
};

struct St_SuperKmer
{
    St_KmerInfo stInfo;
    vector<St_KmerInfo> vSimilar;
    map<uint64_t, int> mpSimilar;
    map<string, St_KmerIndex> mpOrgSimilar;

    //int GetFreqOrg()
    //{
    //    return GetSizeOrg() + stInfo.iNum;
    //}
    void BuildMapOrg(int iGroupIndex);

    /*
    void BuildMap();

    int GetSize();
    */

    int GetSizeOrg();
};

struct St_SimilarKmer
{
    vector<St_KmerInfo> vRegKmer;
    map<string, St_KmerIndex> mpSuperKmer;
};

class  ClsThreadFunc
{
public:
    static void* SimilarDetect(void* args);
};

class ClsMultiThread
{
public:
    ClsMultiThread();

public:
    void FindSimilarKmer(map<string, St_KmerIndex>& mpCombine, vector<St_KmerInfo>& vKmerRegular);
};

#endif // CLSMULTITHREAD_H
