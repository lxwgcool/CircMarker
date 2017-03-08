#include "clsresultcomparison.h"
#include "../../ShareLibrary/clsbasealgorithm.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <unistd.h>

const int IMAXOFFSET = 5; //我们发现有的offset差那么几个bp
//const int CHROMINDEX = 2; // 1 means chromosome 1, 2 means chromosome 2 and etc

ClsResultComparison::ClsResultComparison()
{
    iMyProgramHitNum = 0;
    iCIRIHitNum = 0;
    iFindCircHitNUm = 0;
    iCircExplorerHitNum = 0;

    fMyProgramTP = .0;
    fCIRITP = .0;
    fFindCircTP = .0;
    fCircExplorerTP = .0;
}

//Compare current result with simulation result ===>
void ClsResultComparison::CompareMyResultWithSimulation(string strCurPath, string strSimulationResult,
                                                        int iMaxSupportReads)
{
    cout << "============Compare My_Result With Simulation============" << endl;
    //1: Get Cur Candi
    vector<St_Candidate> vCurCandi;
    ParseBrief(strCurPath, vCurCandi);

    //我们在这里尝试将大于10的结果给过滤掉，看看最后的结果如何 ==>
    for(vector<St_Candidate>::iterator itr = vCurCandi.end() - 1; itr >= vCurCandi.begin(); itr--)
    {
        if(itr->iSupportNum > iMaxSupportReads)
            vCurCandi.erase(itr);
    }
    //<==

    cout << "Current Candiate Sum: " << IntToStr(vCurCandi.size()) << endl;

    //2: Get simulation Candi --> We call simulation caididate as the standard candidate!
    vector<St_Candidate> vStdCandi;
    ParseSimulationResult(strSimulationResult, vStdCandi);
    cout << "Std Candi Sum       : " << IntToStr(vStdCandi.size()) << endl;

    //3: Make Comparison
    //vector<St_Candidate>::iterator itrCur; --> Record the minimum/maximum hit number
    int iMin = 999999;
    int iMax = 0;
    for(vector<St_Candidate>::iterator itr = vStdCandi.end() - 1; itr >= vStdCandi.begin(); itr--)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator itrCur = vCurCandi.begin();
            itrCur != vCurCandi.end(); itrCur++)
        {
            if(CompareTwoCandi(&(*itrCur), &(*itr)))
            {
                //找到了
                itr->iSupportNum = itrCur->iSupportNum;
                //Get minium and maximum hit
                if(iMin > itr->iSupportNum)
                    iMin = itr->iSupportNum;
                if(iMax < itr->iSupportNum)
                    iMax = itr->iSupportNum;
                bFind = true;
                break;
            }
        }

        if(!bFind)
        {
            vStdCandi.erase(itr);
        }

        /*
        itrCur = find(vCurCandi.begin(), vCurCandi.end(), *itr);
        if(itrCur != vCurCandi.end())
        {
            //找到了
            itr->iSupportNum = itrCur->iSupportNum;
            //Get minium and maximum hit
            if(iMin > itr->iSupportNum)
                iMin = itr->iSupportNum;
            if(iMax < itr->iSupportNum)
                iMax = itr->iSupportNum;
        }
        else
        {
            //没找到
            vStdCandi.erase(itr);
        }*/
    }

    //4: Output Result
    cout << "Match Number: " << vStdCandi.size() << endl;
    cout << "TP: " << GetRatio((float)vStdCandi.size() / vCurCandi.size()) << endl;
    cout << "Min Hit Num : " << IntToStr(iMin) << endl;
    cout << "Max Hit Num : " << IntToStr(iMax) << endl;

    cout << "=================================" << endl;
}

void ClsResultComparison::CompareCiriWithSimulation(string strCIRIResult, string strSimulationResult)
{
    cout << "========= Compare Ciri With Simulation =========" << endl;
    //Prase CIRI Result
    vector<St_Candidate> vCiriCandi;
    ParseCIRIResult(strCIRIResult, vCiriCandi);
    cout << "Ciri Candiate Sum: " << IntToStr(vCiriCandi.size()) << endl;

    //Parse simulation result
    vector<St_Candidate> vStdCandi;
    ParseSimulationResult(strSimulationResult, vStdCandi);
    cout << "Std Candi Sum       : " << IntToStr(vStdCandi.size()) << endl;

    //3: Make Comparison
    int iMin = 999999;
    int iMax = 0;
    for(vector<St_Candidate>::iterator itr = vStdCandi.end() - 1; itr >= vStdCandi.begin(); itr--)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator itrCiri = vCiriCandi.begin();
            itrCiri != vCiriCandi.end(); itrCiri++)
        {
            if(CompareTwoCandi(&(*itrCiri), &(*itr)))
            {
                //找到了
                itr->iSupportNum = itrCiri->iSupportNum;
                //Get minium and maximum hit
                if(iMin > itr->iSupportNum)
                    iMin = itr->iSupportNum;
                if(iMax < itr->iSupportNum)
                    iMax = itr->iSupportNum;
                bFind = true;
                break;
            }
        }

        if(!bFind)
        {
            vStdCandi.erase(itr);
        }
    }

    //4: Output Result
    cout << "Match Number: " << vStdCandi.size() << endl;
    cout << "TP: " << GetRatio((float)vStdCandi.size() / vCiriCandi.size()) << endl;
    cout << "Min Hit Num : " << IntToStr(iMin) << endl;
    cout << "Max Hit Num : " << IntToStr(iMax) << endl;

    cout << "=================================" << endl;

}

void ClsResultComparison::CompareFindCircWithSimulation(string strFindCircResult, string strSimulationResult)
{
    cout << "========= Compare Find_Circ With Simulation =========" << endl;
    //Prase CIRI Result
    vector<St_Candidate> vFindCircCandi;
    ParseFindCircResult(strFindCircResult, vFindCircCandi);
    cout << "Find_Circ Candiate Sum: " << IntToStr(vFindCircCandi.size()) << endl;

    //Parse simulation result
    vector<St_Candidate> vStdCandi;
    ParseSimulationResult(strSimulationResult, vStdCandi);
    cout << "Std Candi Sum       : " << IntToStr(vStdCandi.size()) << endl;

    //3: Make Comparison
    int iMin = 999999;
    int iMax = 0;
    for(vector<St_Candidate>::iterator itr = vStdCandi.end() - 1; itr >= vStdCandi.begin(); itr--)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator itrFindCirc = vFindCircCandi.begin();
            itrFindCirc != vFindCircCandi.end(); itrFindCirc++)
        {
            if(CompareTwoCandi(&(*itrFindCirc), &(*itr)))
            {
                //找到了
                itr->iSupportNum = itrFindCirc->iSupportNum;
                //Get minium and maximum hit
                if(iMin > itr->iSupportNum)
                    iMin = itr->iSupportNum;
                if(iMax < itr->iSupportNum)
                    iMax = itr->iSupportNum;
                bFind = true;
                break;
            }
        }

        if(!bFind)
        {
            vStdCandi.erase(itr);
        }
    }

    //4: Output Result
    cout << "Match Number: " << vStdCandi.size() << endl;
    cout << "TP: " << GetRatio((float)vStdCandi.size() / vFindCircCandi.size()) << endl;
    cout << "Min Hit Num : " << IntToStr(iMin) << endl;
    cout << "Max Hit Num : " << IntToStr(iMax) << endl;

    cout << "=================================" << endl;

}

void ClsResultComparison::CompareCircExplorerWithSimulation(string strCIRCexplorerResult, string strSimulationResult)
{
    cout << "========= Compare Circ_Explorer With Simulation =========" << endl;

    //Prase Circ_Explorer Result
    vector<St_Candidate> vCircExpCandi;
    ParseCircExplorerResult(strCIRCexplorerResult, vCircExpCandi);
    cout << "Circ_Explorer Candiate Sum: " << IntToStr(vCircExpCandi.size()) << endl;

    //Parse simulation result
    vector<St_Candidate> vStdCandi;
    ParseSimulationResult(strSimulationResult, vStdCandi);
    cout << "Std Candi Sum       : " << IntToStr(vStdCandi.size()) << endl;

    //3: Make Comparison
    int iMin = 999999;
    int iMax = 0;
    for(vector<St_Candidate>::iterator itr = vStdCandi.end() - 1; itr >= vStdCandi.begin(); itr--)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator itrCircExp = vCircExpCandi.begin();
            itrCircExp != vCircExpCandi.end(); itrCircExp++)
        {
            if(CompareTwoCandi(&(*itrCircExp), &(*itr)))
            {
                //找到了
                itr->iSupportNum = itrCircExp->iSupportNum;
                //Get minium and maximum hit
                if(iMin > itr->iSupportNum)
                    iMin = itr->iSupportNum;
                if(iMax < itr->iSupportNum)
                    iMax = itr->iSupportNum;
                bFind = true;
                break;
            }
        }

        if(!bFind)
        {
            vStdCandi.erase(itr);
        }
    }

    //4: Output Result
    cout << "Match Number: " << vStdCandi.size() << endl;
    cout << "TP: " << GetRatio((float)vStdCandi.size() / vCircExpCandi.size()) << endl;
    cout << "Min Hit Num : " << IntToStr(iMin) << endl;
    cout << "Max Hit Num : " << IntToStr(iMax) << endl;

    cout << "=================================" << endl;
}
//<=====

//=========> Check intersection
void ClsResultComparison::CheckIntersectionOfMyProgram(string strCurResult,
                                                       string strRibominusCur, int iMaxSupportReads)
{    
    cout << endl << "==============Check Intersection of MyProgram==============" << endl;
    if(access(strCurResult.c_str(), 0) != 0 || access(strRibominusCur.c_str(), 0) != 0)
    {
        cout << "One of file Do not existed" << endl;
        return;
    }
    //Parse LC
    //1: Get LC Candi
    vector<St_Candidate> vCurCandi;
    ParseBrief(strCurResult, vCurCandi);        
    cout << "Circular Only based sample: " << IntToStr(vCurCandi.size()) << endl;
    ///Make a filter (filter all of reads which larger than imax hit number)
    for(vector<St_Candidate>::iterator itr = vCurCandi.end() - 1; itr >= vCurCandi.begin(); itr--)
    {
        if(itr->iSupportNum > iMaxSupportReads)
            vCurCandi.erase(itr);
    }
    cout << "Circular Only based sample: " << IntToStr(vCurCandi.size()) << endl;

    //2: Get COnly Candi
    vector<St_Candidate> vRibominusCurCandi;
    ParseBrief(strRibominusCur, vRibominusCurCandi);
    cout << "Linear & Circular based sample: " << IntToStr(vRibominusCurCandi.size()) << endl;
    ofstream ofs;
    ofs.open("./LCInterMyProgram.txt");

    //3: Check intesection --> Output the intersection part
    int iIntersectionNum = 0;
    for(vector<St_Candidate>::iterator itr = vCurCandi.begin(); itr != vCurCandi.end(); itr++)
    {
        for(vector<St_Candidate>::iterator subItr = vRibominusCurCandi.begin();
            subItr != vRibominusCurCandi.end(); subItr++)
        {
            if(CompareTwoCandi(&(*itr), &(*subItr)))
            {
                //output intersection info -->
                ofs << IntToStr(itr->ucChromIndex + 1) << " "
                    << IntToStr(itr->iStartPos) << " "
                    << IntToStr(itr->iEndPos) << " "
                    << "1" << " "  //this is support Num: IntToStr(itr->iSupportNum)
                    << itr->strTag << " "
                    << itr->GetTypeString() << " "
                    << endl;
                //<--

                iIntersectionNum++;
                break;
            }
        }
    }
    ofs.close();
    cout << "Intersection Num: " << IntToStr(iIntersectionNum) << endl;
    cout << "Circular Only based sample ratio    : " << GetRatio((float)iIntersectionNum/vCurCandi.size()) << endl;
    cout << "Linear & Circular based sample ratio: " << GetRatio((float)iIntersectionNum/vRibominusCurCandi.size()) << endl;
}

void ClsResultComparison::CheckIntersectionOfCIRI(string strCIRIResult, string strRibominusCIRI)
{
    cout << endl << "==============Check Intersection of CIRI==============" << endl;
    if(access(strCIRIResult.c_str(), 0) != 0 || access(strRibominusCIRI.c_str(), 0) != 0)
    {
        cout << "One of file Do not existed" << endl;
        return;
    }
    //Parse LC
    //1: Get LC Candi
    vector<St_Candidate> vCIRICandi;
    ParseCIRIResult(strCIRIResult, vCIRICandi);
    cout << "Circular Only based sample: " << IntToStr(vCIRICandi.size()) << endl;

    //2: Get COnly Candi
    vector<St_Candidate> vRibominusCIRICandi;
    ParseCIRIResult(strRibominusCIRI, vRibominusCIRICandi);
    cout << "Linear & Circular based sample: " << IntToStr(vRibominusCIRICandi.size()) << endl;

    ofstream ofs;
    ofs.open("./LCInterCIRI.txt");

    //3: Check intesection
    int iIntersectionNum = 0;
    for(vector<St_Candidate>::iterator itr = vCIRICandi.begin(); itr != vCIRICandi.end(); itr++)
    {
        for(vector<St_Candidate>::iterator subItr = vRibominusCIRICandi.begin();
            subItr != vRibominusCIRICandi.end(); subItr++)
        {
            if(CompareTwoCandi(&(*itr), &(*subItr)))
            {
                //output intersection info -->
                ofs << IntToStr(itr->ucChromIndex + 1) << " "
                    << IntToStr(itr->iStartPos) << " "
                    << IntToStr(itr->iEndPos) << " "
                    << "1" << " "  //this is support Num: IntToStr(itr->iSupportNum)
                    << itr->strTag << " "
                    << itr->GetTypeString() << " "
                    << endl;
                //<--
                iIntersectionNum++;
                break;
            }
        }
    }
    ofs.close();
    cout << "Intersection Num: " << IntToStr(iIntersectionNum) << endl;
    cout << "Circular Only based sample ratio    : " << GetRatio((float)iIntersectionNum/vCIRICandi.size()) << endl;
    cout << "Linear & Circular based sample ratio: " << GetRatio((float)iIntersectionNum/vRibominusCIRICandi.size()) << endl;
}

void ClsResultComparison::CheckIntersectionOfFindCirc(string strFindCircResult, string strRibominusFindCirc)
{
    cout << endl << "==============Check Intersection of Find_Circ==============" << endl;
    if(access(strFindCircResult.c_str(), 0) != 0 || access(strRibominusFindCirc.c_str(), 0) != 0)
    {
        cout << "One of file Do not existed" << endl;
        return;
    }
    //Parse LC
    //1: Get LC Candi
    vector<St_Candidate> vFindCircCandi;
    ParseFindCircResult(strFindCircResult, vFindCircCandi);
    cout << "Circular Only based sample: " << IntToStr(vFindCircCandi.size()) << endl;

    //2: Get COnly Candi
    vector<St_Candidate> vRibominusFindCircCandi;
    ParseFindCircResult(strRibominusFindCirc, vRibominusFindCircCandi);
    cout << "Linear & Circular based sample: " << IntToStr(vRibominusFindCircCandi.size()) << endl;

    ofstream ofs;
    ofs.open("./LCInterFindCirc.txt");

    //3: Check intesection
    int iIntersectionNum = 0;
    for(vector<St_Candidate>::iterator itr = vFindCircCandi.begin(); itr != vFindCircCandi.end(); itr++)
    {
        for(vector<St_Candidate>::iterator subItr = vRibominusFindCircCandi.begin();
            subItr != vRibominusFindCircCandi.end(); subItr++)
        {
            if(CompareTwoCandi(&(*itr), &(*subItr)))
            {
                //output intersection info -->
                ofs << IntToStr(itr->ucChromIndex + 1) << " "
                    << IntToStr(itr->iStartPos) << " "
                    << IntToStr(itr->iEndPos) << " "
                    << "1" << " "  //this is support Num: IntToStr(itr->iSupportNum)
                    << itr->strTag << " "
                    << itr->GetTypeString() << " "
                    << endl;
                //<--
                iIntersectionNum++;
                break;
            }
        }
    }
    ofs.close();
    cout << "Intersection Num: " << IntToStr(iIntersectionNum) << endl;
    cout << "Circular Only based sample ratio    : " << GetRatio((float)iIntersectionNum/vFindCircCandi.size()) << endl;
    cout << "Linear & Circular based sample ratio: "
         << GetRatio((float)iIntersectionNum/vRibominusFindCircCandi.size()) << endl;
}

void ClsResultComparison::CheckIntersectionOfCircExplorer(string strCIRCexplorerResult, string strRibominusCIRCexplorer)
{
    cout << endl << "==============Check Intersection of CIRC_Explorer==============" << endl;
    if(access(strCIRCexplorerResult.c_str(), 0) != 0 || access(strRibominusCIRCexplorer.c_str(), 0) != 0)
    {
        cout << "One of file Do not existed" << endl;
        return;
    }
    //Parse LC
    //1: Get LC Candi
    vector<St_Candidate> vCIRCexplorerCandi;
    ParseCircExplorerResult(strCIRCexplorerResult, vCIRCexplorerCandi);
    cout << "Circular Only based sample: " << IntToStr(vCIRCexplorerCandi.size()) << endl;

    //2: Get COnly Candi
    vector<St_Candidate> vRibominusCIRCexplorerCandi;
    ParseCircExplorerResult(strRibominusCIRCexplorer, vRibominusCIRCexplorerCandi);
    cout << "Linear & Circular based sample: " << IntToStr(vRibominusCIRCexplorerCandi.size()) << endl;

    ofstream ofs;
    ofs.open("./LCInterCIRCexplorer.txt");

    //3: Check intesection
    int iIntersectionNum = 0;
    for(vector<St_Candidate>::iterator itr = vCIRCexplorerCandi.begin(); itr != vCIRCexplorerCandi.end(); itr++)
    {
        for(vector<St_Candidate>::iterator subItr = vRibominusCIRCexplorerCandi.begin();
            subItr != vRibominusCIRCexplorerCandi.end(); subItr++)
        {
            if(CompareTwoCandi(&(*itr), &(*subItr)))
            {
                //output intersection info -->
                ofs << IntToStr(itr->ucChromIndex + 1) << " "
                    << IntToStr(itr->iStartPos) << " "
                    << IntToStr(itr->iEndPos) << " "
                    << "1" << " "  //this is support Num: IntToStr(itr->iSupportNum)
                    << itr->strTag << " "
                    << itr->GetTypeString() << " "
                    << endl;
                //<--
                iIntersectionNum++;
                break;
            }
        }
    }
    ofs.close();
    cout << "Intersection Num: " << IntToStr(iIntersectionNum) << endl;
    cout << "Circular Only based sample ratio    : " << GetRatio((float)iIntersectionNum/vCIRCexplorerCandi.size()) << endl;
    cout << "Linear & Circular based sample ratio: "
         << GetRatio((float)iIntersectionNum/vRibominusCIRCexplorerCandi.size()) << endl;
}

void ClsResultComparison::GetMostSupportedCandiAndCheckIntersection(string strCurResult,
                                                                    string strCIRIResult,
                                                                    string strFindCircResult,
                                                                    string strCIRCexplorerResult)
{
    cout << endl << "=============GetMostSupportedCandiAndCheckIntersection===========" << endl;
    //1: Get the most supported Candidate
    /// 1: Parse the intersection result of "strCurResult"
    vector<St_Candidate> vCurResult;
    ParseBrief("./LCInterMyProgram.txt", vCurResult);
    EraseDuplicate(vCurResult);

    /// 2: Parse the intersection result of "strCIRIResult"
    vector<St_Candidate> vCIRI;
    ParseBrief("./LCInterCIRI.txt", vCIRI);
    EraseDuplicate(vCIRI);

    /// 3: Parse the intersection result of "strFindCircResult"
    vector<St_Candidate> vFindCirc;
    ParseBrief("./LCInterFindCirc.txt", vFindCirc);
    EraseDuplicate(vFindCirc);

    /// 4: Parse the intersection result of "strCIRCexplorerResult"
    vector<St_Candidate> vCIRCexplorer;
    ParseBrief("./LCInterCIRCexplorer.txt", vCIRCexplorer);
    EraseDuplicate(vCIRCexplorer);

    /// 5: Collect most reliable candidate
    vector<St_Candidate> vReliableCandi;
    //////(1) Consider My_Program
    for(vector<St_Candidate>::iterator itr = vCurResult.begin(); itr != vCurResult.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vReliableCandi.begin();
            subItr != vReliableCandi.end(); subItr++)
        {
            if(CompareTwoCandi(&(*itr), &(*subItr)))
            {
                subItr->iSupportNum++;
                bFind = true;
                break;
            }
        }
        if(!bFind) // 没找到，是新的
        {
            vReliableCandi.push_back(*itr);
        }
    }

    //////(2) Consider CIRI
    for(vector<St_Candidate>::iterator itr = vCIRI.begin(); itr != vCIRI.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vReliableCandi.begin();
            subItr != vReliableCandi.end(); subItr++)
        {
            if(CompareTwoCandi(&(*itr), &(*subItr)))
            {
                subItr->iSupportNum++;
                bFind = true;
                break;
            }
        }
        if(!bFind) // 没找到，是新的
        {
            vReliableCandi.push_back(*itr);
        }
    }

    //////(3) Consider FindCirc
    for(vector<St_Candidate>::iterator itr = vFindCirc.begin(); itr != vFindCirc.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vReliableCandi.begin();
            subItr != vReliableCandi.end(); subItr++)
        {
            if(CompareTwoCandi(&(*itr), &(*subItr)))
            {
                subItr->iSupportNum++;
                bFind = true;
                break;
            }
        }
        if(!bFind) // 没找到，是新的
        {
            vReliableCandi.push_back(*itr);
        }
    }

    //////(4) Consider CIRCexplorer
    for(vector<St_Candidate>::iterator itr = vCIRCexplorer.begin(); itr != vCIRCexplorer.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vReliableCandi.begin();
            subItr != vReliableCandi.end(); subItr++)
        {
            if(CompareTwoCandi(&(*itr), &(*subItr)))
            {
                subItr->iSupportNum++;
                bFind = true;
                break;
            }
        }
        if(!bFind) // 没找到，是新的
        {
            vReliableCandi.push_back(*itr);
        }
    }

    ///////(5)Only Select the Reiable Candidate with support num >=2
    for(vector<St_Candidate>::iterator itr = vReliableCandi.end() - 1; itr >= vReliableCandi.begin(); itr--)
    {
        if(itr->iSupportNum < 2)
        {
            vReliableCandi.erase(itr);
        }
    }
    cout << "vReliableCandi Size (support >=2): " << IntToStr(vReliableCandi.size()) << endl;

    //2: Check Interesection
    ///(1) Parse strCurResult
    vector<St_Candidate> vCurCandi;
    ParseBrief(strCurResult, vCurCandi);
    EraseDuplicate(vCurCandi);
    int iCurSupport = 0;
    for(vector<St_Candidate>::iterator itr = vReliableCandi.begin(); itr != vReliableCandi.end(); itr++)
    {
        for(vector<St_Candidate>::iterator subItr = vCurCandi.begin(); subItr != vCurCandi.end(); subItr++)
        {
            if(CompareTwoCandi(&(*itr), &(*subItr)))
            {
                iCurSupport++;
                break;
            }
        }
    }
    cout << "My_Program and Reliable Full Set Comparison: "
         << IntToStr(iCurSupport) << " " << GetRatio(float(iCurSupport)/vReliableCandi.size()) << endl;

    ///(2) Parse strCIRIResult
    vector<St_Candidate> vCIRICandi;
    ParseCIRIResult(strCIRIResult, vCIRICandi);
    EraseDuplicate(vCIRICandi);
    int iCIRISupport = 0;
    for(vector<St_Candidate>::iterator itr = vReliableCandi.begin(); itr != vReliableCandi.end(); itr++)
    {
        for(vector<St_Candidate>::iterator subItr = vCIRICandi.begin(); subItr != vCIRICandi.end(); subItr++)
        {
            if(CompareTwoCandi(&(*itr), &(*subItr)))
            {
                iCIRISupport++;
                break;
            }
        }
    }
    cout << "CIRI and Reliable Full Set Comparison: "
         << IntToStr(iCIRISupport) << " " << GetRatio(float(iCIRISupport)/vReliableCandi.size()) << endl;

    ///(3) Parse strFindCircResult
    vector<St_Candidate> vFindCircCandi;
    ParseFindCircResult(strFindCircResult, vFindCircCandi);
    EraseDuplicate(vFindCircCandi);
    int iFindCircSupport = 0;
    for(vector<St_Candidate>::iterator itr = vReliableCandi.begin(); itr != vReliableCandi.end(); itr++)
    {
        for(vector<St_Candidate>::iterator subItr = vFindCircCandi.begin();
            subItr != vFindCircCandi.end(); subItr++)
        {
            if(CompareTwoCandi(&(*itr), &(*subItr)))
            {
                iFindCircSupport++;
                break;
            }
        }
    }
    cout << "FindCirc and Reliable Full Set Comparison: "
         << IntToStr(iFindCircSupport) << " " << GetRatio(float(iFindCircSupport)/vReliableCandi.size()) << endl;

    ///(4) Parse strCIRCexplorerResult
    vector<St_Candidate> vCIRCexplorerCandi;
    ParseCircExplorerResult(strCIRCexplorerResult, vCIRCexplorerCandi);
    EraseDuplicate(vCIRCexplorerCandi);
    int iCIRCexplorerSupport = 0;
    for(vector<St_Candidate>::iterator itr = vReliableCandi.begin(); itr != vReliableCandi.end(); itr++)
    {
        for(vector<St_Candidate>::iterator subItr = vCIRCexplorerCandi.begin();
            subItr != vCIRCexplorerCandi.end(); subItr++)
        {
            if(CompareTwoCandi(&(*itr), &(*subItr)))
            {
                iCIRCexplorerSupport++;
                break;
            }
        }
    }
    cout << "CIRCexplorer and Reliable Full Set Comparison: "
         << IntToStr(iCIRCexplorerSupport) << " " << GetRatio(float(iCIRCexplorerSupport)/vReliableCandi.size()) << endl;
}

void ClsResultComparison::EraseDuplicate(vector<St_Candidate>& vResult)
{
    if(vResult.empty())
        return;

    for(vector<St_Candidate>::iterator itr = vResult.end() - 1; itr >= vResult.begin(); itr--)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vResult.begin(); subItr <= itr-1; subItr++)
        {
            if(*itr == *subItr)
            {
                bFind = true;
                break;
            }
        }
        if(bFind)
        {
            vResult.erase(itr);
        }
    }
}

//<=========

void ClsResultComparison::CompareGTFWithCircBaseResult(vector<St_Row_Chrom>& vChrom,
                                                       string strCircBasePath)
{
    //Get the information of CircBase Result
    vector<St_Candidate> vCandi;
    ParseCircBaseResult(strCircBasePath, vCandi);

    //record the unhit std candidate
    string strUnHitCandiFile = "./un_hit_candi.txt";
    string strHitCandiFile = "./hit_candi.txt";
    string strChr1ValidCandi = "./chr1_valid_candi.txt";
    ofstream ofsUnHit;
    ofstream ofsHit;
    ofstream ofsChr1Candi;
    ofsUnHit.open(strUnHitCandiFile.c_str());
    ofsHit.open(strHitCandiFile.c_str());
    ofsChr1Candi.open(strChr1ValidCandi.c_str());

    //Make Comparison
    int iHitNum = 0;
    int iValidCandi = 0;
    int iChr1ValidCandi = 0;
    int iChr1TotalCandi = 0;    

    for(vector<St_Candidate>::iterator itr = vCandi.begin(); itr != vCandi.end(); itr++)
    {
        bool bHitBothStartAndEnd = false;
        int iCheckResult = CheckBothHitExonBoundary(vChrom, &(*itr));

        if(iCheckResult == 0) // 这个candidate所属的chromosome根本在gtf中没有找到
            continue;

        if(iCheckResult == 1) //两边都比对上了
        {
            bHitBothStartAndEnd = true;
        }
        iValidCandi++;

        if(itr->ucChromIndex == 0)
            iChr1TotalCandi++;

        if(bHitBothStartAndEnd)
        {
            iHitNum++;
            ofsHit << "chr" << IntToStr(itr->ucChromIndex + 1) << "\t"
                   << "start: " << IntToStr(itr->iStartPos) << "\t"
                   << "end: " << IntToStr(itr->iEndPos) << "\t"
                   << "Direction: " << (itr->bRC ? "-" : "+") << endl;

            if(itr->ucChromIndex == 0) // this is the first chromosone in reference
            {
                iChr1ValidCandi++;
                ofsChr1Candi << "chr" << IntToStr(itr->ucChromIndex + 1) << "\t"
                             << "start: " << IntToStr(itr->iStartPos) << "\t"
                             << "end: " << IntToStr(itr->iEndPos) << "\t"
                             << "Direction: " << (itr->bRC ? "-" : "+") << endl;
            }
        }
        else
        {
            ofsUnHit << "chr" << IntToStr(itr->ucChromIndex + 1) << "\t"
                << "start: " << IntToStr(itr->iStartPos) << "\t"
                << "end: " << IntToStr(itr->iEndPos) << "\t"
                << "Direction: " << (itr->bRC ? "-" : "+") << endl;
        }
    }

    ofsHit.close();
    ofsUnHit.close();
    ofsChr1Candi.close();

    //Output Result
    cout << "===========CompareGTFWithCircBaseResult===========" << endl;
    cout << "CircBase Valid Candidate Sum  : " << IntToStr(iValidCandi) << endl;
    cout << "Hit Bundary with GTF Sum      : " << IntToStr(iHitNum) << endl;
    cout << "Hit Chromosome 1 Candidata    : " << IntToStr(iChr1ValidCandi) << endl;
    cout << "Hit Chr 1 Total Candidata     : " << IntToStr(iChr1TotalCandi) << endl << endl;
}

void ClsResultComparison::CompareMyResultWithCircBaseResult(vector<St_Row_Chrom>& vChrom,
                                                            string strCircBasePath,
                                                            string strCurResultPath,
                                                            int iCurChromIndex)
{
    cout << "===========CompareGTFWithCircBaseResult===========" << endl;
    if(access(strCircBasePath.c_str(), 0) != 0 || access(strCurResultPath.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    ///Step1: Get my Result
    vector<St_Candidate> vCurResult;
    ParseBrief(strCurResultPath, vCurResult);
    cout << "Brief Cur Num: " << IntToStr(vCurResult.size()) << endl;

    for(vector<St_Candidate>::iterator itr = vCurResult.begin(); itr < vCurResult.end(); itr++) //this logic is correct
    {
        for(vector<St_Candidate>::iterator subItr = vCurResult.end() - 1; subItr > itr; subItr--)
        {
            if(itr->ucChromIndex == subItr->ucChromIndex) // 同一个 chromosome
            {
                if( (itr->iStartPos == subItr->iStartPos && itr->iEndPos == subItr->iEndPos) ||
                    (itr->iEndPos == subItr->iStartPos && itr->iStartPos == subItr->iEndPos))
                {
                    itr->iSupportNum += subItr->iSupportNum;
                    vCurResult.erase(subItr);
                }
            }
        }
    }
    cout << "Final Cur Result Num: " << IntToStr(vCurResult.size()) << endl;

    ///Step2: Get the valid candi in CircBase: (need to combine both vChrom and CircBase)
    vector<St_Candidate> vCandi;
    ParseCircBaseResult(strCircBasePath, vCandi);
    //我们在这里暂时只看chromosome 1
    ofstream ofsChr1Candi;
    ofsChr1Candi.open("./chr1_valid_std_candi.txt");
    vector<St_Candidate> vCandiForSpeciChrom;
    int iChr1StdCandiNum = 0;
    for(vector<St_Candidate>::iterator itr = vCandi.begin(); itr != vCandi.end(); itr++)
    {
        if(itr->ucChromIndex != iCurChromIndex - 1) //我们在这里只考虑chromosome1的standard candidate
            continue;

        iChr1StdCandiNum++;
        int iCheckResult = CheckBothHitExonBoundary(vChrom, &(*itr));

        if(iCheckResult == 1) // both map
        {
            vCandiForSpeciChrom.push_back(*itr);

            ofsChr1Candi << "chr" << IntToStr(itr->ucChromIndex + 1) << "\t"
                         << "start: " << IntToStr(itr->iStartPos) << "\t"
                         << "end: " << IntToStr(itr->iEndPos) << "\t"
                         << "Direction: " << (itr->bRC ? "-" : "+") << endl;

        }
    }

    ofsChr1Candi.close();

    cout << "Total Chr1 Std Candi Number: " << IntToStr(iChr1StdCandiNum) << endl;
    cout << "Valid Chr1 Std Candi Number: " << IntToStr(vCandiForSpeciChrom.size()) << endl << endl;

    ///Step3: Make the comparison: here "start pos" is smaller than "end pos"
    ofstream ofs;
    ofs.open("./HittedCurCandi.txt");
    ofstream ofsUnHit;
    ofsUnHit.open("./Un_Hit_Cur_Candi.txt");
    ofstream ofsUnHitStdCandi;
    ofsUnHitStdCandi.open("./Un_Hit_Std_Candi.txt");

    int iHitNum = 0;
    vector<St_Candidate> vHitCandi;
    //-->
    int iSelfCircNum = 0;
    int iRegCircNum = 0;
    int iSelfHit = 0;
    int iRegHit = 0;

    int iEE = 0;
    int iSS = 0;
    int iSE = 0;
    int iES = 0;
    int iLG = 0;

    int iEESum = 0;
    int iSSSum = 0;
    int iSESum = 0;
    int iESSum = 0;
    int iLGSum = 0;
    //<--
    for(vector<St_Candidate>::iterator itr = vCurResult.begin(); itr != vCurResult.end(); itr++)
    {
        if(itr->enCircType == ctSelf)
            iSelfCircNum++;
        else if(itr->enCircType == ctRegular)
        {
            iRegCircNum++;                

            if(itr->strTag == "EE")
                iEESum++;
            else if(itr->strTag == "SS")
                iSSSum++;
            else if(itr->strTag == "SE")
                iSESum++;
            else if(itr->strTag == "ES")
                iESSum++;
            else if(itr->strTag == "LG")
                iLGSum++;
            else
            {}
        }

        bool bFind = false;      
        for(vector<St_Candidate>::iterator subItr = vCandiForSpeciChrom.begin();
            subItr != vCandiForSpeciChrom.end(); subItr++)
        {
            bFind = CompareTwoCandi(&(*itr), &(*subItr));

            if(bFind)
            {                                
                vHitCandi.push_back(*subItr);
                break;
            }
        }

        if(bFind)
        {
            //-->
            if(itr->enCircType == ctSelf)
                iSelfHit++;
            else if(itr->enCircType == ctRegular)
            {
                iRegHit++;
                if(itr->strTag == "EE")
                    iEE++;
                else if(itr->strTag == "SS")
                    iSS++;
                else if(itr->strTag == "SE")
                    iSE++;
                else if(itr->strTag == "ES")
                    iES++;
                else if(itr->strTag == "LG")
                    iLG++;
                else
                {}
            }
            //<--

            iHitNum++;            
            //Save the hitted cur candi
            ofs << IntToStr(itr->ucChromIndex + 1) << " "
                << IntToStr(itr->iStartPos) << " "
                << IntToStr(itr->iEndPos) << " "
                << IntToStr(itr->iSupportNum) << " "
                << itr->strTag << " "
                << itr->GetTypeString() << " "
                << endl;
        }
        else
        {
            ofsUnHit << "chr" << IntToStr(itr->ucChromIndex + 1) << "\t"
                     << IntToStr(itr->iStartPos) << "\t"
                     << IntToStr(itr->iEndPos) << "\t"
                     << IntToStr(itr->iSupportNum) << "\t"
                     << itr->strTag << "\t"
                     << itr->GetTypeString() << "\t"
                     << endl;
        }
    }

    iMyProgramHitNum = iHitNum;
    fMyProgramTP = (float)iHitNum/vCurResult.size();
    cout << "Total Cur Candi : " << IntToStr(vCurResult.size()) << endl;
    cout << "Hitted Cur Candi: " << IntToStr(iHitNum) << "\t"
         << GetRatio((float)iHitNum/vCurResult.size()) << endl;

    cout << "Self Circ   : " << IntToStr(iSelfCircNum) << " --> Hit: " << IntToStr(iSelfHit) << endl;
    cout << "Regular Circ: " << IntToStr(iRegCircNum) << " --> Hit: " << IntToStr(iRegHit) << endl << endl;
    cout << "EE: " << IntToStr(iEESum) << " -->" << IntToStr(iEE) << endl;
    cout << "SS: " << IntToStr(iSSSum) << " -->" << IntToStr(iSS) << endl;
    cout << "SE: " << IntToStr(iSESum) << " -->" << IntToStr(iSE) << endl;
    cout << "ES: " << IntToStr(iESSum) << " -->" << IntToStr(iES) << endl;
    cout << "LG: " << IntToStr(iLGSum) << " -->" << IntToStr(iLG) << endl;

    //Save the unhit std candidate
    for(vector<St_Candidate>::iterator itr = vCandiForSpeciChrom.begin();
        itr != vCandiForSpeciChrom.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vHitCandi.begin();
            subItr != vHitCandi.end(); subItr++)
        {
            if(*itr == *subItr)
            {
                bFind = true;
                break;
            }
        }
        if(!bFind)
        {
            ofsUnHitStdCandi << "chr" << IntToStr(itr->ucChromIndex + 1) << "\t"
                             << IntToStr(itr->iStartPos) << "\t"
                             << IntToStr(itr->iEndPos) << "\t"
                             << IntToStr(itr->iSupportNum) << endl;
        }
    }

    ofs.close();
    ofsUnHit.close();
    ofsUnHitStdCandi.close();
}

void ClsResultComparison::CompareCiriWithCircBaseResult(vector<St_Row_Chrom>& vChrom,
                                                        string strCircBasePath,
                                                        string strCiriResultPath,
                                                        int iCurChromIndex)
{
    if(access(strCircBasePath.c_str(), 0) != 0 || access(strCiriResultPath.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    //1: Parse CircBase Result
    vector<St_Candidate> vCandiCircBase;
    ParseCircBaseResult(strCircBasePath, vCandiCircBase);

    //2: Just get the chromsome one's candidate for CircBase
    //我们在这里暂时只看chromosome 1
    vector<St_Candidate> vCandiForSpeciChrom;
    for(vector<St_Candidate>::iterator itr = vCandiCircBase.begin(); itr != vCandiCircBase.end(); itr++)
    {
        if(itr->ucChromIndex != iCurChromIndex - 1) //我们在这里只考虑chromosome1的standard candidate
            continue;
        int iCheckResult = CheckBothHitExonBoundary(vChrom, &(*itr));

        if(iCheckResult == 1) // both map
        {
            vCandiForSpeciChrom.push_back(*itr);
        }
    }

    //3: Parse CiriResult
    vector<St_Candidate> vCandiCIRI;
    ParseCIRIResult(strCiriResultPath, vCandiCIRI);

    //4: Make Comparison
    ofstream ofs;
    ofs.open("./CIRI_HittedCandi.txt");
    int iHitNum = 0;
    for(vector<St_Candidate>::iterator itr = vCandiCIRI.begin(); itr != vCandiCIRI.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vCandiForSpeciChrom.begin();
            subItr != vCandiForSpeciChrom.end(); subItr++)
        {
            bFind = CompareTwoCandi(&(*itr), &(*subItr));

            if(bFind)
            {
                ofs << IntToStr(itr->ucChromIndex) << " "
                    << IntToStr(itr->iStartPos) << " "
                    << IntToStr(itr->iEndPos) << " "
                    << IntToStr(itr->iSupportNum) << " "
                    << itr->strTag << " "
                    << itr->GetTypeString() << " "
                    << endl;
                iHitNum++;
                break;
            }
        }
    }

    iCIRIHitNum = iHitNum;
    fCIRITP = (float)iHitNum / vCandiCIRI.size();

    cout << "------------------CompareCiriWithCircBaseResult------------------" << endl;
    cout << "vCandiForSpeciChrom Size: " << IntToStr(vCandiForSpeciChrom.size()) << endl;
    cout << "vCandiCIRI Size    : " << IntToStr(vCandiCIRI.size()) << endl;
    cout << "Hit Num            : " << IntToStr(iHitNum) << "\t"
         << GetRatio((float)iHitNum / vCandiCIRI.size()) << endl;
    ofs.close();
}

void ClsResultComparison::CompareFindCircWithCircBaseResult( vector<St_Row_Chrom>& vChrom,
                                                             string strCircBasePath,
                                                             string strFindCircPath,
                                                             int iCurChromIndex)
{
    if(access(strCircBasePath.c_str(), 0) != 0 || access(strFindCircPath.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    //1: Parse CircBase Result
    vector<St_Candidate> vCandiCircBase;
    ParseCircBaseResult(strCircBasePath, vCandiCircBase);

    ///Just get the chromsome one's candidate for CircBase
    //我们在这里暂时只看chromosome 1
    vector<St_Candidate> vCandiForSpeciChrom;
    for(vector<St_Candidate>::iterator itr = vCandiCircBase.begin(); itr != vCandiCircBase.end(); itr++)
    {
        if(itr->ucChromIndex != iCurChromIndex - 1) //我们在这里只考虑chromosome1的standard candidate
            continue;
        //int iCheckResult = CheckBothHitExonBoundary(vChrom, &(*itr));

        //if(iCheckResult == 1) // both map
        {
            vCandiForSpeciChrom.push_back(*itr);
        }
    }

    //2: Parse Find_Circ Result
    vector<St_Candidate> vCandiFindCirc;
    ParseFindCircResult(strFindCircPath, vCandiFindCirc);  // do not -1

    /*
    ///--->Additional Code
    ///Erase the candidate which do not constructed with the boundary of exon
    for(vector<St_Candidate>::iterator itr = vCandiFindCirc.end()-1; itr >= vCandiFindCirc.begin(); itr--)
    {
        if(itr->ucChromIndex != iCurChromIndex - 1) //我们在这里只考虑chromosome1的standard candidate
            continue;
        int iCheckResult = CheckBothHitExonBoundary(vChrom, &(*itr));

        if(iCheckResult != 1) // both map
        {
            vCandiFindCirc.erase(itr);
        }
    }
    ///<----
    /// */

    //3: Make comparison
    ofstream ofs;
    ofs.open("./FindCirc_HittedCandi.txt");
    int iHitNum = 0;
    for(vector<St_Candidate>::iterator itr = vCandiFindCirc.begin(); itr != vCandiFindCirc.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vCandiForSpeciChrom.begin();
            subItr != vCandiForSpeciChrom.end(); subItr++)
        {
            bFind = CompareTwoCandi(&(*itr), &(*subItr));

            if(bFind)
            {
                ofs << IntToStr(itr->ucChromIndex) << " "
                    << IntToStr(itr->iStartPos) << " "
                    << IntToStr(itr->iEndPos) << " "
                    << IntToStr(itr->iSupportNum) << " "
                    << itr->strTag << " "
                    << itr->GetTypeString() << " "
                    << endl;
                iHitNum++;
                break;
            }
        }
    }

    iFindCircHitNUm = iHitNum;
    fFindCircTP = (float)iHitNum / vCandiFindCirc.size();

    cout << "------------------CompareFindCircWithCircBaseResult------------------" << endl;
    cout << "vCandiForSpeciChrom Size: " << IntToStr(vCandiForSpeciChrom.size()) << endl;
    cout << "vCandiFindCirc Size     : " << IntToStr(vCandiFindCirc.size()) << endl;
    cout << "Hit Num                 : " << IntToStr(iHitNum) << "\t"
         << GetRatio((float)iHitNum / vCandiFindCirc.size()) << endl;

    ofs.close();
}

void ClsResultComparison::CompareCircExplorerWithCircBaseResult(vector<St_Row_Chrom>& vChrom,
                                                                string strCircBasePath,
                                                                string strCircExplorerPath,
                                                                int iCurChromIndex)
{
    if(access(strCircBasePath.c_str(), 0) != 0 || access(strCircExplorerPath.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    //1: Parse CircBase Result
    vector<St_Candidate> vCandiCircBase;
    ParseCircBaseResult(strCircBasePath, vCandiCircBase);

    ///Just get the chromsome one's candidate for CircBase
    //我们在这里暂时只看chromosome 1
    vector<St_Candidate> vCandiForSpeciChrom;
    for(vector<St_Candidate>::iterator itr = vCandiCircBase.begin(); itr != vCandiCircBase.end(); itr++)
    {
        if(itr->ucChromIndex != iCurChromIndex - 1) //我们在这里只考虑chromosome1的standard candidate
            continue;
        int iCheckResult = CheckBothHitExonBoundary(vChrom, &(*itr));

        if(iCheckResult == 1) // both map
        {
            vCandiForSpeciChrom.push_back(*itr);
        }
    }

    //2: Parse Find_Circ Result
    vector<St_Candidate> vCandiCircExplorer;
    ParseCircExplorerResult(strCircExplorerPath, vCandiCircExplorer);

    //3: Make comparison
    ofstream ofs;
    ofs.open("./CircExplorer_HittedCandi.txt");
    int iHitNum = 0;
    for(vector<St_Candidate>::iterator itr = vCandiCircExplorer.begin();
        itr != vCandiCircExplorer.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vCandiForSpeciChrom.begin();
            subItr != vCandiForSpeciChrom.end(); subItr++)
        {
            bFind = CompareTwoCandi(&(*itr), &(*subItr));

            if(bFind)
            {
                ofs << IntToStr(itr->ucChromIndex) << " "
                    << IntToStr(itr->iStartPos) << " "
                    << IntToStr(itr->iEndPos) << " "
                    << IntToStr(itr->iSupportNum) << " "
                    << itr->strTag << " "
                    << itr->GetTypeString() << " "
                    << endl;

                iHitNum++;
                break;
            }
        }
    }

    iCircExplorerHitNum = iHitNum;
    fCircExplorerTP = (float)iHitNum / vCandiCircExplorer.size();

    cout << "------------------CompareCIRCexplorerWithCircBaseResult------------------" << endl;
    cout << "vCandiForSpeciChrom Size: " << IntToStr(vCandiForSpeciChrom.size()) << endl;
    cout << "vCandiCircExplorer Size : " << IntToStr(vCandiCircExplorer.size()) << endl;
    cout << "Hit Num                 : " << IntToStr(iHitNum) << "\t"
         << GetRatio((float)iHitNum / vCandiCircExplorer.size()) << endl;

    ofs.close();
}

void ClsResultComparison::SummaryCompareResult()
{
    cout << endl;
    //Get the Norm Hit Num;
    int iSum = iMyProgramHitNum + iCIRIHitNum + iFindCircHitNUm + iCircExplorerHitNum;
    cout << "Hit_Num Avg Ratio My Program, CIRI, FindCirc, CircExplorer: "
         << GetRatio((float)iMyProgramHitNum/iSum) << ", "
         << GetRatio((float)iCIRIHitNum/iSum) << ", "
         << GetRatio((float)iFindCircHitNUm/iSum) << ", "
         << GetRatio((float)iCircExplorerHitNum/iSum) << endl;

    float fSum = fMyProgramTP + fCIRITP + fFindCircTP + fCircExplorerTP;
    cout << "TP Avg Ratio My Program, CIRI, FindCirc, CircExplorer: "
         << GetRatio((float)fMyProgramTP/fSum) << ", "
         << GetRatio((float)fCIRITP/fSum) << ", "
         << GetRatio((float)fFindCircTP/fSum) << ", "
         << GetRatio((float)fCircExplorerTP/fSum) << endl;

    cout << "Avg Ratio Sum My Program, CIRI, FindCirc, CircExplorer: "
         << FloatToStr((float)iMyProgramHitNum/iSum + (float)fMyProgramTP/fSum, 3) << ", "
         << FloatToStr((float)iCIRIHitNum/iSum + (float)fCIRITP/fSum, 3) << ", "
         << FloatToStr((float)iFindCircHitNUm/iSum + (float)fFindCircTP/fSum, 3) << ", "
         << FloatToStr((float)iCircExplorerHitNum/iSum + (float)fCircExplorerTP/fSum, 3) << endl;

    cout << endl;
    //Get the Avg Hit Num;
    int iMax = iMyProgramHitNum > iCIRIHitNum ? iMyProgramHitNum : iCIRIHitNum;
    iMax = iMax > iFindCircHitNUm ? iMax : iFindCircHitNUm;
    iMax = iMax > iCircExplorerHitNum ? iMax : iCircExplorerHitNum;
    cout << "Hit_Num Norm Ratio My Program, CIRI, FindCirc, CircExplorer: "
         << GetRatio((float)iMyProgramHitNum/iMax) << ", "
         << GetRatio((float)iCIRIHitNum/iMax) << ", "
         << GetRatio((float)iFindCircHitNUm/iMax) << ", "
         << GetRatio((float)iCircExplorerHitNum/iMax) << endl;

    float fMax = fMyProgramTP > fCIRITP ? fMyProgramTP : fCIRITP;
    fMax = fMax > fFindCircTP ? fMax : fFindCircTP;
    fMax = fMax > fCircExplorerTP ? fMax : fCircExplorerTP;
    cout << "TP Norm Ratio My Program, CIRI, FindCirc, CircExplorer: "
         << GetRatio((float)fMyProgramTP/fMax) << ", "
         << GetRatio((float)fCIRITP/fMax) << ", "
         << GetRatio((float)fFindCircTP/fMax) << ", "
         << GetRatio((float)fCircExplorerTP/fMax) << endl;

    cout << "Norm Ratio Sum My Program, CIRI, FindCirc, CircExplorer: "
         << FloatToStr((float)iMyProgramHitNum/iMax + (float)fMyProgramTP/fMax, 3) << ", "
         << FloatToStr((float)iCIRIHitNum/iMax + (float)fCIRITP/fMax, 3) << ", "
         << FloatToStr((float)iFindCircHitNUm/iMax + (float)fFindCircTP/fMax, 3) << ", "
         << FloatToStr((float)iCircExplorerHitNum/iMax + (float)fCircExplorerTP/fMax, 3) << endl;

    cout << endl;
}

//*************************************************************************
//************For My_Program Compare with Other three**********************
//*************************************************************************
void ClsResultComparison::CheckIntersetBTCIRIAndMyProgram(string strCIRI, string strMyProgram)
{
    if(access(strCIRI.c_str(), 0) != 0 || access(strMyProgram.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    //1: Parse CIRI Result
    vector<St_Candidate> vCandiCIRI;
    ParseBrief(strCIRI, vCandiCIRI);

    //2: Parse My Program Result
    vector<St_Candidate> vCandiMyProgram;
    ParseBrief(strMyProgram, vCandiMyProgram);

    ofstream ofs;
    ofs.open("Only_CIRI_MP.txt");
    int iOnlyNum = 0;
    for(vector<St_Candidate>::iterator itr = vCandiCIRI.begin();
        itr != vCandiCIRI.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vCandiMyProgram.begin();
            subItr != vCandiMyProgram.end(); subItr++)
        {
            bFind = CompareTwoCandi(&(*itr), &(*subItr));
            if(bFind)
            {
                break;
            }
        }
        if(!bFind)
        {
            ofs << "chr" << IntToStr(itr->ucChromIndex) << "\t"
                << IntToStr(itr->iStartPos) << "\t"
                << IntToStr(itr->iEndPos) << "\t"
                << IntToStr(itr->iSupportNum) << endl;
            iOnlyNum++;
        }
    }
    ofs.close();
    cout << "Only_CIRI_MP Num: " << IntToStr(iOnlyNum)
         << " --> " << GetRatio(1-(float)iOnlyNum/vCandiCIRI.size())<<endl;
}

void ClsResultComparison::CheckIntersetBTFindCircAndMyProgram(string strFindCirc, string strMyProgram)
{
    if(access(strFindCirc.c_str(), 0) != 0 || access(strMyProgram.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    //1: Parse FindCirc Result
    vector<St_Candidate> vCandiFindCirc;
    ParseBrief(strFindCirc, vCandiFindCirc);  // 这里需要的是hitted 的find circ的路径

    //2: Parse My Program Result
    vector<St_Candidate> vMyProgram;
    ParseBrief(strMyProgram, vMyProgram);

    ofstream ofs;
    ofs.open("Only_FindCirc_MP.txt");
    int iOnlyNum = 0;
    for(vector<St_Candidate>::iterator itr = vCandiFindCirc.begin();
        itr != vCandiFindCirc.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vMyProgram.begin();
            subItr != vMyProgram.end(); subItr++)
        {
            bFind = CompareTwoCandi(&(*itr), &(*subItr));
            if(bFind)
            {
                break;
            }
        }
        if(!bFind)
        {
            ofs << "chr" << IntToStr(itr->ucChromIndex) << "\t"
                << IntToStr(itr->iStartPos) << "\t"
                << IntToStr(itr->iEndPos) << "\t"
                << IntToStr(itr->iSupportNum) << endl;
            iOnlyNum++;
        }
    }
    ofs.close();
    cout << "Only_FindCirc_MP Num: " << IntToStr(iOnlyNum)
         << " --> " << GetRatio(1-(float)iOnlyNum/vCandiFindCirc.size())<<endl;
}

void ClsResultComparison::CheckIntersetBTCircExplorerAndMyProgram(string strCircExplorer, string strMyProgram)
{
    if(access(strCircExplorer.c_str(), 0) != 0 || access(strMyProgram.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    //1: Parse FindCirc Result
    vector<St_Candidate> vCandiCircExplorer;
    ParseBrief(strCircExplorer, vCandiCircExplorer);  // 这里需要的是hitted 的find circ的路径

    //2: Parse My Program Result
    vector<St_Candidate> vMyProgram;
    ParseBrief(strMyProgram, vMyProgram);

    ofstream ofs;
    ofs.open("Only_CircExplorer_MP.txt");
    int iOnlyNum = 0;
    for(vector<St_Candidate>::iterator itr = vCandiCircExplorer.begin();
        itr != vCandiCircExplorer.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vMyProgram.begin();
            subItr != vMyProgram.end(); subItr++)
        {
            bFind = CompareTwoCandi(&(*itr), &(*subItr));
            if(bFind)
            {
                break;
            }
        }
        if(!bFind)
        {
            ofs << "chr" << IntToStr(itr->ucChromIndex) << "\t"
                << IntToStr(itr->iStartPos) << "\t"
                << IntToStr(itr->iEndPos) << "\t"
                << IntToStr(itr->iSupportNum) << endl;
            iOnlyNum++;
        }
    }
    ofs.close();
    cout << "Only_CircExplorer_MP Num: " << IntToStr(iOnlyNum)
         << " --> " << GetRatio(1-(float)iOnlyNum/vCandiCircExplorer.size())<<endl;
}


//*************************************************************************
//************For CIRI Compare with Other three**********************
//*************************************************************************
void ClsResultComparison::CheckIntersetBTMyProgramAndCIRI(string strMyProgram, string strCIRI)
{
    if(access(strMyProgram.c_str(), 0) != 0 || access(strCIRI.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    //1: Parse My Program Result
    vector<St_Candidate> vCandiMyProgram;
    ParseBrief(strMyProgram, vCandiMyProgram);

    //2: Parse CIRI Result
    vector<St_Candidate> vCandiCIRI;
    ParseBrief(strCIRI, vCandiCIRI);

    ofstream ofs;
    ofs.open("Only_MyProgram_CIRI.txt");
    int iOnlyNum = 0;
    for(vector<St_Candidate>::iterator itr = vCandiMyProgram.begin();
        itr != vCandiMyProgram.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vCandiCIRI.begin();
            subItr != vCandiCIRI.end(); subItr++)
        {
            bFind = CompareTwoCandi(&(*itr), &(*subItr));
            if(bFind)
            {
                break;
            }
        }
        if(!bFind)
        {
            ofs << "chr" << IntToStr(itr->ucChromIndex) << "\t"
                << IntToStr(itr->iStartPos) << "\t"
                << IntToStr(itr->iEndPos) << "\t"
                << IntToStr(itr->iSupportNum) << endl;
            iOnlyNum++;
        }
    }
    ofs.close();
    cout << "Only_MyProgram_CIRI Num: " << IntToStr(iOnlyNum)
         << " --> " << GetRatio(1-(float)iOnlyNum/vCandiMyProgram.size())<<endl;
}

void ClsResultComparison::CheckIntersetBTFindCircAndCIRI(string strFindCirc, string strCIRI)
{
    if(access(strFindCirc.c_str(), 0) != 0 || access(strCIRI.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    //1: Parse Find_Circ Result
    vector<St_Candidate> vCandiFindCirc;
    ParseBrief(strFindCirc, vCandiFindCirc);

    //2: Parse CIRI Result
    vector<St_Candidate> vCandiCIRI;
    ParseBrief(strCIRI, vCandiCIRI);

    ofstream ofs;
    ofs.open("Only_FindCirc_CIRI.txt");
    int iOnlyNum = 0;
    for(vector<St_Candidate>::iterator itr = vCandiFindCirc.begin();
        itr != vCandiFindCirc.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vCandiCIRI.begin();
            subItr != vCandiCIRI.end(); subItr++)
        {
            bFind = CompareTwoCandi(&(*itr), &(*subItr));
            if(bFind)
            {
                break;
            }
        }
        if(!bFind)
        {
            ofs << "chr" << IntToStr(itr->ucChromIndex) << "\t"
                << IntToStr(itr->iStartPos) << "\t"
                << IntToStr(itr->iEndPos) << "\t"
                << IntToStr(itr->iSupportNum) << endl;
            iOnlyNum++;
        }
    }
    ofs.close();
    cout << "Only_FindCirc_CIRI Num: " << IntToStr(iOnlyNum)
         << " --> " << GetRatio(1-(float)iOnlyNum/vCandiFindCirc.size())<<endl;
}

void ClsResultComparison::CheckIntersetBTCircExplorerAndCIRI(string strCircExplorer, string strCIRI)
{
    if(access(strCircExplorer.c_str(), 0) != 0 || access(strCIRI.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    //1: Parse Circ_Explorer Result
    vector<St_Candidate> vCandiCircExplorer;
    ParseBrief(strCircExplorer, vCandiCircExplorer);

    //2: Parse CIRI Result
    vector<St_Candidate> vCandiCIRI;
    ParseBrief(strCIRI, vCandiCIRI);

    ofstream ofs;
    ofs.open("Only_CircExplorer_CIRI.txt");
    int iOnlyNum = 0;
    for(vector<St_Candidate>::iterator itr = vCandiCircExplorer.begin();
        itr != vCandiCircExplorer.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vCandiCIRI.begin();
            subItr != vCandiCIRI.end(); subItr++)
        {
            bFind = CompareTwoCandi(&(*itr), &(*subItr));
            if(bFind)
            {
                break;
            }
        }
        if(!bFind)
        {
            ofs << "chr" << IntToStr(itr->ucChromIndex) << "\t"
                << IntToStr(itr->iStartPos) << "\t"
                << IntToStr(itr->iEndPos) << "\t"
                << IntToStr(itr->iSupportNum) << endl;
            iOnlyNum++;
        }
    }
    ofs.close();
    cout << "Only_CircExplorer_CIRI Num: " << IntToStr(iOnlyNum)
         << " --> " << GetRatio(1-(float)iOnlyNum/vCandiCircExplorer.size())<<endl;
}

//*************************************************************************
//************For Find_Circ Compare with Other three**********************
//*************************************************************************
void ClsResultComparison::CheckIntersetBTMyProgramAndFindCirc(string strMyProgram, string strFindCirc)
{
    if(access(strMyProgram.c_str(), 0) != 0 || access(strFindCirc.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    //1: Parse My Program Result
    vector<St_Candidate> vCandiMyProgram;
    ParseBrief(strMyProgram, vCandiMyProgram);

    //2: Parse Find_Circ Result
    vector<St_Candidate> vCandiFindCirc;
    ParseBrief(strFindCirc, vCandiFindCirc);

    ofstream ofs;
    ofs.open("Only_MyProgram_FC.txt");
    int iOnlyNum = 0;
    for(vector<St_Candidate>::iterator itr = vCandiMyProgram.begin();
        itr != vCandiMyProgram.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vCandiFindCirc.begin();
            subItr != vCandiFindCirc.end(); subItr++)
        {
            bFind = CompareTwoCandi(&(*itr), &(*subItr));
            if(bFind)
            {
                break;
            }
        }
        if(!bFind)
        {
            ofs << "chr" << IntToStr(itr->ucChromIndex) << "\t"
                << IntToStr(itr->iStartPos) << "\t"
                << IntToStr(itr->iEndPos) << "\t"
                << IntToStr(itr->iSupportNum) << endl;
            iOnlyNum++;
        }
    }
    ofs.close();
    cout << "Only_MyProgram_FC Num: " << IntToStr(iOnlyNum)
         << " --> " << GetRatio(1-(float)iOnlyNum/vCandiMyProgram.size())<<endl;
}

void ClsResultComparison::CheckIntersetBTCIRIAndFindCirc(string strCIRI, string strFindCirc)
{
    if(access(strCIRI.c_str(), 0) != 0 || access(strFindCirc.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    //1: Parse CIRI Result
    vector<St_Candidate> vCandiCIRI;
    ParseBrief(strCIRI, vCandiCIRI);

    //2: Parse Find_Circ Result
    vector<St_Candidate> vCandiFindCirc;
    ParseBrief(strFindCirc, vCandiFindCirc);

    ofstream ofs;
    ofs.open("Only_CIRI_FC.txt");
    int iOnlyNum = 0;
    for(vector<St_Candidate>::iterator itr = vCandiCIRI.begin();
        itr != vCandiCIRI.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vCandiFindCirc.begin();
            subItr != vCandiFindCirc.end(); subItr++)
        {
            bFind = CompareTwoCandi(&(*itr), &(*subItr));
            if(bFind)
            {
                break;
            }
        }
        if(!bFind)
        {
            ofs << "chr" << IntToStr(itr->ucChromIndex) << "\t"
                << IntToStr(itr->iStartPos) << "\t"
                << IntToStr(itr->iEndPos) << "\t"
                << IntToStr(itr->iSupportNum) << endl;
            iOnlyNum++;
        }
    }
    ofs.close();
    cout << "Only_CIRI_FC Num: " << IntToStr(iOnlyNum)
         << " --> " << GetRatio(1-(float)iOnlyNum/vCandiCIRI.size())<<endl;
}

void ClsResultComparison::CheckIntersetBTCircExplorerAndFindCirc(string strCircExplorer, string strFindCirc)
{
    if(access(strCircExplorer.c_str(), 0) != 0 || access(strFindCirc.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    //1: Parse Circ_Explorer Result
    vector<St_Candidate> vCandiCircExplorer;
    ParseBrief(strCircExplorer, vCandiCircExplorer);

    //2: Parse Find_Circ Result
    vector<St_Candidate> vCandiFindCirc;
    ParseBrief(strFindCirc, vCandiFindCirc);

    ofstream ofs;
    ofs.open("Only_CircExplorer_FC.txt");
    int iOnlyNum = 0;
    for(vector<St_Candidate>::iterator itr = vCandiCircExplorer.begin();
        itr != vCandiCircExplorer.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vCandiFindCirc.begin();
            subItr != vCandiFindCirc.end(); subItr++)
        {
            bFind = CompareTwoCandi(&(*itr), &(*subItr));
            if(bFind)
            {
                break;
            }
        }
        if(!bFind)
        {
            ofs << "chr" << IntToStr(itr->ucChromIndex) << "\t"
                << IntToStr(itr->iStartPos) << "\t"
                << IntToStr(itr->iEndPos) << "\t"
                << IntToStr(itr->iSupportNum) << endl;
            iOnlyNum++;
        }
    }
    ofs.close();
    cout << "Only_CircExplorer_FC Num: " << IntToStr(iOnlyNum)
         << " --> " << GetRatio(1-(float)iOnlyNum/vCandiCircExplorer.size())<<endl;
}


//*************************************************************************
//************For Circ_Explorer Compare with Other three**********************
//*************************************************************************
void ClsResultComparison::CheckIntersetBTMyProgramAndCircExplorer(string strMyProgram, string strCircExplorer)
{
    if(access(strMyProgram.c_str(), 0) != 0 || access(strCircExplorer.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    //1: Parse My Program Result
    vector<St_Candidate> vCandiMyProgram;
    ParseBrief(strMyProgram, vCandiMyProgram);

    //2: Parse Circ_Explorer Result
    vector<St_Candidate> vCandiCircExplorer;
    ParseBrief(strCircExplorer, vCandiCircExplorer);

    ofstream ofs;
    ofs.open("Only_MyProgram_CE.txt");
    int iOnlyNum = 0;
    for(vector<St_Candidate>::iterator itr = vCandiMyProgram.begin();
        itr != vCandiMyProgram.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vCandiCircExplorer.begin();
            subItr != vCandiCircExplorer.end(); subItr++)
        {
            bFind = CompareTwoCandi(&(*itr), &(*subItr));
            if(bFind)
            {
                break;
            }
        }
        if(!bFind)
        {
            ofs << "chr" << IntToStr(itr->ucChromIndex) << "\t"
                << IntToStr(itr->iStartPos) << "\t"
                << IntToStr(itr->iEndPos) << "\t"
                << IntToStr(itr->iSupportNum) << endl;
            iOnlyNum++;
        }
    }
    ofs.close();
    cout << "Only_MyProgram_CE Num: " << IntToStr(iOnlyNum)
         << " --> " << GetRatio(1-(float)iOnlyNum/vCandiMyProgram.size())<<endl;
}

void ClsResultComparison::CheckIntersetBTCIRIAndCircExplorer(string strFindCirc, string strCircExplorer)
{
    if(access(strFindCirc.c_str(), 0) != 0 || access(strCircExplorer.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    //1: Parse Find_Circ Result
    vector<St_Candidate> vCandiFindCirc;
    ParseBrief(strFindCirc, vCandiFindCirc);

    //2: Parse Circ_Explorer Result
    vector<St_Candidate> vCandiCircExplorer;
    ParseBrief(strCircExplorer, vCandiCircExplorer);

    ofstream ofs;
    ofs.open("Only_FindCirc_CE.txt");
    int iOnlyNum = 0;
    for(vector<St_Candidate>::iterator itr = vCandiFindCirc.begin();
        itr != vCandiFindCirc.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vCandiCircExplorer.begin();
            subItr != vCandiCircExplorer.end(); subItr++)
        {
            bFind = CompareTwoCandi(&(*itr), &(*subItr));
            if(bFind)
            {
                break;
            }
        }
        if(!bFind)
        {
            ofs << "chr" << IntToStr(itr->ucChromIndex) << "\t"
                << IntToStr(itr->iStartPos) << "\t"
                << IntToStr(itr->iEndPos) << "\t"
                << IntToStr(itr->iSupportNum) << endl;
            iOnlyNum++;
        }
    }
    ofs.close();
    cout << "Only_FindCirc_CE Num: " << IntToStr(iOnlyNum)
         << " --> " << GetRatio(1-(float)iOnlyNum/vCandiFindCirc.size())<<endl;
}

void ClsResultComparison::CheckIntersetBTFindCircAndCircExplorer(string strCIRI, string strCircExplorer)
{
    if(access(strCIRI.c_str(), 0) != 0 || access(strCircExplorer.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    //1: Parse CIRI Result
    vector<St_Candidate> vCandiCIRI;
    ParseBrief(strCIRI, vCandiCIRI);

    //2: Parse Circ_Explorer Result
    vector<St_Candidate> vCandiCircExplorer;
    ParseBrief(strCircExplorer, vCandiCircExplorer);

    ofstream ofs;
    ofs.open("Only_CIRI_CE.txt");
    int iOnlyNum = 0;
    for(vector<St_Candidate>::iterator itr = vCandiCIRI.begin();
        itr != vCandiCIRI.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vCandiCircExplorer.begin();
            subItr != vCandiCircExplorer.end(); subItr++)
        {
            bFind = CompareTwoCandi(&(*itr), &(*subItr));
            if(bFind)
            {
                break;
            }
        }
        if(!bFind)
        {
            ofs << "chr" << IntToStr(itr->ucChromIndex) << "\t"
                << IntToStr(itr->iStartPos) << "\t"
                << IntToStr(itr->iEndPos) << "\t"
                << IntToStr(itr->iSupportNum) << endl;
            iOnlyNum++;
        }
    }
    ofs.close();
    cout << "Only_CIRI_CE Num: " << IntToStr(iOnlyNum)
         << " --> " << GetRatio(1-(float)iOnlyNum/vCandiCIRI.size())<<endl;
}

///---------------Base Parse Function-----------------
/// --------------------------------------------------
/// --------------------------------------------------
// Brief.txt (same directory as executable file)
void ClsResultComparison::ParseBrief(string strCurPath, vector<St_Candidate>& vCandi)
{
    ifstream ifs;
    ifs.open(strCurPath.c_str());
    vCandi.clear();
    string strLine = "";
    St_Candidate stCandi; //candidate atom

    while(!ifs.eof()) // check if reached the end
    {
        getline(ifs, strLine);
        //1: For Chrom Index
        int iStart = 0;
        int iEnd = strLine.find(" ");
        int iLen = iEnd - iStart;
        stCandi.ucChromIndex = atoi(strLine.substr(iStart, iLen).c_str());

        //2: Start Pos
        iStart = iEnd + 1;
        iEnd = strLine.find(" ", iStart);
        iLen = iEnd - iStart;
        int iFirstValue = atoi(strLine.substr(iStart, iLen).c_str());

        //3: End Pos
        iStart = iEnd + 1;
        iEnd = strLine.find(" ", iStart);
        iLen = iEnd - iStart;
        int iSecondValue = atoi(strLine.substr(iStart, iLen).c_str());

        stCandi.iStartPos = iFirstValue <= iSecondValue ? iFirstValue : iSecondValue;
        stCandi.iEndPos = iFirstValue >= iSecondValue ? iFirstValue : iSecondValue;

        //4:Support Reads Count
        iStart = iEnd + 1;
        iEnd = strLine.find(" ", iStart);
        iLen = iEnd - iStart;
        stCandi.iSupportNum = atoi(strLine.substr(iStart, iLen).c_str());

        //5: Tag
        iStart = iEnd + 1;
        iEnd = strLine.find(" ", iStart);
        iLen = iEnd - iStart;
        stCandi.strTag = strLine.substr(iStart, iLen);

        //6: Circular Type
        iStart = iEnd + 1;
        iEnd = strLine.length();
        iLen = iEnd - iStart;
        string strType = strLine.substr(iStart, iLen);
        if(strType == "S")
            stCandi.enCircType = ctSelf;
        else if(strType == "R")
            stCandi.enCircType = ctRegular;

        //6: Save Current Candidate
        //if(stCandi.iSupportNum < 5) // half coverage
            vCandi.push_back(stCandi);
    }
    ifs.close();
}

// Standard Result  --------->Do it tomorrow
// we only do the testing based on chr1
void ClsResultComparison::ParseSimulationResult(string strStdPath, vector<St_Candidate>& vCandi)
{
    ifstream ifs;
    ifs.open(strStdPath.c_str());
    vCandi.clear();
    string strLine = "";
    St_Candidate stCandi; //candidate atom

    while(!ifs.eof()) // check if reached the end
    {
        getline(ifs, strLine);
        //Get Chromosone Name
        int iStartPos = 0;
        int iEndPos = strLine.find('\t', iStartPos);
        int iLen = iEndPos - iStartPos;
        string strChromName = strLine.substr(iStartPos, iLen);
        if(strChromName == "chr1")
            stCandi.ucChromIndex = 0;
        else
            stCandi.ucChromIndex = 128; // Max Value

        //Get First value
        iEndPos = strLine.rfind('\t');
        iStartPos = strLine.rfind('\t', iEndPos - 1) + 1;
        iLen = iEndPos - iStartPos;
        int iFirstValue = atoi(strLine.substr(iStartPos, iLen).c_str());

        //Get Second Value
        iStartPos = strLine.rfind('\t') + 1;
        iLen = strLine.length() - iStartPos;
        int iSecondValue = atoi(strLine.substr(iStartPos, iLen).c_str());

        stCandi.iStartPos = iFirstValue <= iSecondValue ? iFirstValue : iSecondValue;
        stCandi.iEndPos = iFirstValue >= iSecondValue ? iFirstValue : iSecondValue;

        //在此处将重复的进行排除
        bool bExisted = false;
        for(vector<St_Candidate>::iterator itr = vCandi.begin(); itr != vCandi.end(); itr++)
        {
            if( *itr == stCandi)
            {
                bExisted = true;
                break;
            }
        }

        if(!bExisted) // only record it when it is not existed
            vCandi.push_back(stCandi);
    }
    ifs.close();
}

//CircBase Standard Circular RNA
void ClsResultComparison::ParseCircBaseResult(string strCircBasePath, vector<St_Candidate>& vCandi)
{
    //Parse CircBase Reuslt
    /* File Structure
     * 0: chromosone name
     * 1: start
     * 2: end
     *    (在这里好像都是从小到大的)
     * 3: direction: + or -
     */
    ifstream ifs;
    ifs.open(strCircBasePath.c_str());
    vCandi.clear();
    string strLine = "";
    St_Candidate stCandi; //candidate atom

    while(!ifs.eof()) // check if reached the end
    {
        getline(ifs, strLine);

        if(strLine == "")
            continue;

        if(strLine.length() > 0 && strLine.substr(0, 1) == "#")
            continue;

        int iStartPos = 0;
        int iEndPos = strLine.find('\t', iStartPos);
        int iLen = iEndPos - iStartPos;

        //The 0 --> chromosone name
        string strChrName = strLine.substr(iStartPos, iLen);
        iStartPos = 3;
        iLen = strChrName.length() - 3; //
        //注意，这里我们还是从0开始，这样我们可以直接通过下标去找相应的chromosone里面的exons的信息
        stCandi.ucChromIndex = atoi(strChrName.substr(iStartPos, iLen).c_str()) - 1;

        //The 1 --> Start position
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;
        stCandi.iStartPos = atoi(strLine.substr(iStartPos, iLen).c_str());

        //The 2 --> End position
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;
        stCandi.iEndPos = atoi(strLine.substr(iStartPos, iLen).c_str());

        //The 3 --> Direction
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;
        string strTag = strLine.substr(iStartPos, iLen);
        if(strTag == "-")
            stCandi.bRC = true;
        else
            stCandi.bRC = false;

        //在此处将重复的进行排除
        bool bExisted = false;
        for(vector<St_Candidate>::iterator itr = vCandi.begin(); itr != vCandi.end(); itr++)
        {
            if( *itr == stCandi)
            {
                bExisted = true;
                break;
            }
        }

        if(!bExisted) // only record it when it is not existed
            vCandi.push_back(stCandi);

    }
    ifs.close();
}

// CIRI Result
void ClsResultComparison::ParseCIRIResult(string strCIRIPath, vector<St_Candidate>& vCandi)
{
    vCandi.clear();
    ifstream ifs;
    ifs.open(strCIRIPath.c_str(), ios::in);
    string strLine = "";
    St_Candidate stCandi;
    while(!ifs.eof()) // check if reached the end
    {
        getline(ifs, strLine);
        int iPos = strLine.find('\t');
        if(iPos < 0)
            continue;

        string strValue = strLine.substr(0, iPos);
        int iPosFirstSplit = strValue.find(':');
        int iPosSecondSplit = strValue.find('|');
        if(iPosFirstSplit < 0 || iPosSecondSplit < 0)
            continue;

        //Get Real value
        stCandi.ucChromIndex = atoi(strValue.substr(0, iPosFirstSplit).c_str());
        stCandi.iStartPos = atoi(strValue.substr(iPosFirstSplit+1,
                                                  iPosSecondSplit - iPosFirstSplit - 1).c_str());
        stCandi.iEndPos = atoi(strValue.substr(iPosSecondSplit + 1, strValue.length() -
                                                                iPosSecondSplit - 1).c_str());
        if(strLine.find('+') >= 0)
            stCandi.bRC = false;
        else
            stCandi.bRC = true;

        vCandi.push_back(stCandi);
    }
    ifs.close();
}

//CIRCexplorer Result
void ClsResultComparison::ParseCircExplorerResult(string strCIRCexplorerPath,
                                                  vector<St_Candidate>& vCandi)
{
    vCandi.clear();
    ifstream ifs;
    ifs.open(strCIRCexplorerPath.c_str());
    string strLine = "";
    St_Candidate stCandi;

    while(!ifs.eof())
    {
        getline(ifs, strLine);
        int iTagPos = strLine.find('\t');
        if(iTagPos < 0)
            continue;

        //1: For chromosome
        int iStart = 0;
        int iLen = iTagPos - iStart;
        stCandi.ucChromIndex = atoi(strLine.substr(iStart, iLen).c_str()) - 1; //注意这里我们没有-1

        //2: For Start Pos
        iStart = iTagPos + 1;
        iTagPos = strLine.find('\t', iStart);
        iLen = iTagPos  - iStart;
        stCandi.iStartPos = atoi(strLine.substr(iStart, iLen).c_str());

        //For the second value
        iStart = iTagPos + 1;
        iTagPos = strLine.find('\t', iStart);
        iLen = iTagPos  - iStart;
        stCandi.iEndPos = atoi(strLine.substr(iStart, iLen).c_str());

        //For Direction
        if(strLine.find('+', iStart) >= 0)
            stCandi.bRC = false;
        else
            stCandi.bRC = true;

        vCandi.push_back(stCandi);
    }
    ifs.close();
}

//Find_Circ Result
void ClsResultComparison::ParseFindCircResult(string strFindCircPath, vector<St_Candidate>& vCandi)
{
    vCandi.clear();
    ifstream ifs;
    ifs.open(strFindCircPath.c_str());
    string strLine = "";
    St_Candidate stCandi;
    while(!ifs.eof())
    {
        getline(ifs, strLine);
        int iTagPos = strLine.find('\t');
        if(iTagPos < 0)
            continue;

        //1: For chromosome
        int iStart = 0;
        int iLen = iTagPos - iStart;
        stCandi.ucChromIndex = atoi(strLine.substr(iStart, iLen).c_str()) - 1; //注意这里我们没有-1 -->We need -1 to coincide with others

        //2: For Start Pos
        iStart = iTagPos + 1;
        iTagPos = strLine.find('\t', iStart);
        iLen = iTagPos  - iStart;
        stCandi.iStartPos = atoi(strLine.substr(iStart, iLen).c_str());

        //For the second value
        iStart = iTagPos + 1;
        iTagPos = strLine.find('\t', iStart);
        iLen = iTagPos  - iStart;
        stCandi.iEndPos = atoi(strLine.substr(iStart, iLen).c_str());

        //For Direction
        if(strLine.find('+', iStart) >= 0)
            stCandi.bRC = false;
        else
            stCandi.bRC = true;

        vCandi.push_back(stCandi);
    }
    ifs.close();
}

//Sub-Functions
int ClsResultComparison::CheckBothHitExonBoundary(vector<St_Row_Chrom>& vChrom, St_Candidate* pCandi)
{
    St_Row_Chrom *pRawChrom = NULL;
    for(vector<St_Row_Chrom>::iterator itrChrom = vChrom.begin();
        itrChrom != vChrom.end(); itrChrom++)
    {
        if(pCandi->ucChromIndex + 1 == atoi(itrChrom->strName.c_str()))  // here, we just need to add 1 to coincide the rule of index
        {
            pRawChrom = &(*itrChrom);
            break;
        }
    }

    if(pRawChrom == NULL)
        return 0;

    bool bHitBothStartAndEnd = false;

    for(vector<St_Raw_Gene>::iterator itrRG = pRawChrom->vRG.begin();
        itrRG != pRawChrom->vRG.end(); itrRG++)
    {
        for(vector<St_Raw_Transcript>::iterator itrRT = itrRG->vRT.begin();
            itrRT != itrRG->vRT.end(); itrRT++)
        {
            bool bFindStart = false;
            bool bFindEnd = false;
            for(vector<St_Raw_Exon>::iterator itrExon = itrRT->vRExon.begin();
                itrExon != itrRT->vRExon.end(); itrExon++)
            {
                if(!bFindStart)
                {
                    if(abs(pCandi->iStartPos - itrExon->iStart) <= IMAXOFFSET)
                        bFindStart = true;
                    else if(abs(pCandi->iStartPos - itrExon->iEnd) <= IMAXOFFSET)
                        bFindStart = true;
                }

                if(!bFindEnd)
                {
                    if(abs(pCandi->iEndPos - itrExon->iStart) <= IMAXOFFSET)
                        bFindEnd = true;
                    else if(abs(pCandi->iEndPos - itrExon->iEnd) <= IMAXOFFSET)
                        bFindEnd = true;
                }

                if(bFindStart && bFindEnd)
                {
                    bHitBothStartAndEnd = true;
                    break;
                }
            }
            if(bHitBothStartAndEnd)
                break;
        }
        if(bHitBothStartAndEnd)
            break;
    }

    if(bHitBothStartAndEnd)
        return 1;
    else
        return -1;
}

//General function: compare two candidate
bool ClsResultComparison::CompareTwoCandi(St_Candidate* pCandi1, St_Candidate* pCandi2)
{
    bool bFind = false;
    int iStartCur = -1;
    int iEndCur = -1;
    if(pCandi1->iStartPos <= pCandi1->iEndPos)
    {
        iStartCur = pCandi1->iStartPos;
        iEndCur = pCandi1->iEndPos;
    }
    else
    {
        iStartCur = pCandi1->iEndPos;
        iEndCur = pCandi1->iStartPos;
    }

    int iStartStd = -1;
    int iEndStd = -1;

    if(pCandi2->iStartPos <= pCandi2->iEndPos)
    {
        iStartStd = pCandi2->iStartPos;
        iEndStd = pCandi2->iEndPos;
    }
    else
    {
        iStartStd = pCandi2->iEndPos;
        iEndStd = pCandi2->iStartPos;
    }

    //Make the comparison
    if(abs(iStartCur - iStartStd) <= IMAXOFFSET &&
       abs(iEndCur - iEndStd) <= IMAXOFFSET)
    {
        bFind = true;
    }
    return bFind;
}
