#include <iostream>
#include "../../ShareLibrary/clsgtfparse.h"
#include "clsconfig.h"
#include "../../ShareLibrary/clskmeralgorithm.h"
#include "../../ShareLibrary/clsbasealgorithm.h"
#include "clskmertable.h"
#include "clsfindcandidate.h"
#include "clsresultcomparison.h"
#include "clstroubleshoot.h"
#include "string.h"
#include <unistd.h>
using namespace std;

void Test(string strFaPath);

//#define ONLY_DO_COMPARISON
//#define CUR_COMPARISON
//#define SIMULATION_COMPARISON
//#define CHECK_INTERSECTION

void FindCircRNA(St_Config& stConfig, ClsKmerTable* pKT, ClsFindCandidate* pFindCandi,
                 vector<St_Fastq>& vFastq, vector<St_Row_Chrom>& vChrom, vector<St_Fasta>& vFasta);

int main(int argc, char **argv)
{
    if(argc == 1 ||
       strcmp(argv[1], "-h") == 0 ||
       strcmp(argv[1], "--help") == 0)
    {
        cout << "**********************************" << endl;
        cout << "How to use CircMarker" << endl;
        cout << "./CircRnaDetectDraft ./config.ini" << endl;
        cout << "**********************************" << endl;
        return 0;
    }

    //1: Read Ini
    ClsConfig* pConfig = new ClsConfig();
    St_Config stConfig;
    pConfig->ReadConfig(stConfig, argv[1]);
    bool bOk = true;
    if(!pConfig->CheckConfig(stConfig))
    {
        bOk = false;
    }
    delete pConfig;
    pConfig = NULL;

    //Check if File exists
    //1: Check ref
    if(::access(stConfig.strRefPath.c_str(), 0) != 0)
    {
        bOk = false;
        cout << "Error: Reference does not exist!" << endl;
    }

    //2: Check fastq
    if(stConfig.strReads1Path != "" &&
       ::access(stConfig.strReads1Path.c_str(), 0) != 0)
    {
        bOk = false;
        cout << "Error: Reads1 does not exist!" << endl;
    }

    if(stConfig.strReads2Path != "" &&
       ::access(stConfig.strReads2Path.c_str(), 0) != 0)
    {
        bOk = false;
        cout << "Error: Reads2 does not exist!" << endl;
    }

    //3: Check gtf
    if(::access(stConfig.strGtfPath.c_str(), 0) != 0)
    {
        bOk = false;
        cout << "Error: Annotation file (GTF) does not exist!" << endl;
    }

    if(!bOk)
    {
        return 0;
    }

    //Test(stConfig.strRefPath);

    //2: Read GTF---------------------------
    //   Get the information: Chromosn, Gene, Transcript, and Exon
    vector<St_Row_Chrom> vChrom;
    ClsGTFParse* pGTFParse = new ClsGTFParse();
    pGTFParse->Init(stConfig.strGtfPath, stConfig.strRefPath, KMERLEN,
                    stConfig.iReadsLen, stConfig.fKmerRatio);
    pGTFParse->ReadGTF(vChrom); //Done
//#ifndef ONLY_DO_COMPARISON
    pGTFParse->GetTagValue(vChrom);
//#endif
    delete pGTFParse;
    pGTFParse = NULL;

    cout << "The first vChrom Size: " << vChrom.size() << endl;

    //3.0 Read Reference File
    //Read Fasta:
    ClsFastaReader* pFastaReader = new ClsFastaReader();
    vector<St_Fasta> vFasta;
    pFastaReader->ReadFastaRegular(stConfig.strRefPath, vFasta);
    delete pFastaReader;
    pFastaReader = NULL;

    //Check if reference naming rule (chromosome) is the same as gtf naming rule
    if(vChrom.empty() || vFasta.empty())
    {
        cout << "Error: Invalid Annotation or Reference File (chromosome list is empty)!" << endl;
        return 0;
    }
    else
    {
        string strGtfChromExample = vChrom[0].strName;
        bool bFind = false;
        for(vector<St_Fasta>::iterator itrRef = vFasta.begin(); itrRef != vFasta.end(); itrRef++)
        {
            string strRefName = "";
            if(itrRef->strName.find(' ') == string::npos)
                strRefName = itrRef->strName;
            else
                strRefName = itrRef->strName.substr(0, itrRef->strName.find(' '));

            if(strRefName == strGtfChromExample)
            {
                bFind = true;
                break;
            }
        }
        if(!bFind)
        {
            cout << "Error: Ther naming rule of chromosome in annotation file and reference is different (please use NCBI/EMSEMBL sequence level style)!" << endl;
            return 0;
        }
    }

    //3: Create Kmer Table-----------------
    ClsKmerTable* pKT = new ClsKmerTable();    

    //4: Check Hitting
    ClsFindCandidate* pFindCandi = new ClsFindCandidate();

    //1: Read Reads
    vector<St_Fastq> vFastq;
    pFindCandi->AssembleReads(stConfig.strReads1Path, stConfig.strReads2Path, vFastq);


#ifndef ONLY_DO_COMPARISON
    FindCircRNA(stConfig, pKT, pFindCandi, vFastq, vChrom, vFasta);

#endif
//#ifdef ONLY_DO_COMPARISON
    //5: simple testing --> for pre-verification --> Go!!
    ClsResultComparison* pResultComparison = new ClsResultComparison();

#ifdef SIMULATION_COMPARISON
    pResultComparison->CompareMyResultWithSimulation(stConfig.strCurResult,
                                                     stConfig.strSimulationResult, stConfig.iMaxSupportReads);
    pResultComparison->CompareCiriWithSimulation(stConfig.strCIRIResult, stConfig.strSimulationResult);
    pResultComparison->CompareFindCircWithSimulation(stConfig.strFindCircResult, stConfig.strSimulationResult);
    pResultComparison->CompareCircExplorerWithSimulation(stConfig.strCIRCexplorerResult, stConfig.strSimulationResult);
#endif


#ifdef CHECK_INTERSECTION
    pResultComparison->CheckIntersectionOfMyProgram(stConfig.strCurResult,
                                                    stConfig.strRibominusCur,
                                                    stConfig.iMaxSupportReads); // simulation result will be used to save the Circular Only based sample
    pResultComparison->CheckIntersectionOfCIRI(stConfig.strCIRIResult, stConfig.strRibominusCIRI); //  ...
    pResultComparison->CheckIntersectionOfFindCirc(stConfig.strFindCircResult, stConfig.strRibominusFindCirc); //  ...
    pResultComparison->CheckIntersectionOfCircExplorer(stConfig.strCIRCexplorerResult, stConfig.strRibominusCIRCexplorer); //  ...

    pResultComparison->GetMostSupportedCandiAndCheckIntersection(stConfig.strCurResult,
                                                                 stConfig.strCIRIResult,
                                                                 stConfig.strFindCircResult,
                                                                 stConfig.strCIRCexplorerResult);
#endif

    //pResultComparison->CompareGTFWithCircBaseResult(vChrom, stConfig.strCircBaseResult);

#ifdef CUR_COMPARISON
    pResultComparison->CompareMyResultWithCircBaseResult(vChrom,
                                                         stConfig.strCircBaseResult,
                                                         stConfig.strCurResult,
                                                         stConfig.iCurChromIndex);

    pResultComparison->CompareCiriWithCircBaseResult(vChrom, stConfig.strCircBaseResult,
                                                     stConfig.strCIRIResult,
                                                     stConfig.iCurChromIndex);

    pResultComparison->CompareFindCircWithCircBaseResult(vChrom, stConfig.strCircBaseResult,
                                                         stConfig.strFindCircResult,
                                                         stConfig.iCurChromIndex);

    pResultComparison->CompareCircExplorerWithCircBaseResult(vChrom, stConfig.strCircBaseResult,
                                                             stConfig.strCIRCexplorerResult,
                                                             stConfig.iCurChromIndex);
    pResultComparison->SummaryCompareResult();
#endif

#ifdef CUR_COMPARISON    
    cout << "************ My Result Compare with CIRC, Find_Circ, and Circ_Explorer ************" << endl;
    cout << "************ ********************************************************* ************" << endl;
    pResultComparison->CheckIntersetBTCIRIAndMyProgram(stConfig.strCIRIHit, stConfig.strCurHit);


    pResultComparison->CheckIntersetBTFindCircAndMyProgram(stConfig.strFindCircHit, stConfig.strCurHit);
    pResultComparison->CheckIntersetBTCircExplorerAndMyProgram(stConfig.strCircExplorerHit, stConfig.strCurHit);

    // 比对 CIRI VS My program
    //             Find_Circ
    //             Circ_Explorer
    cout << endl << endl;
    cout << "************ CIRI Compare with My Result, Find_Circ, and Circ_Explorer ************" << endl;
    cout << "************ ********************************************************* ************" << endl;
    pResultComparison->CheckIntersetBTMyProgramAndCIRI(stConfig.strCurHit, stConfig.strCIRIHit);
    pResultComparison->CheckIntersetBTFindCircAndCIRI(stConfig.strFindCircHit, stConfig.strCIRIHit);
    pResultComparison->CheckIntersetBTCircExplorerAndCIRI(stConfig.strCircExplorerHit, stConfig.strCIRIHit);

    // 比对 Find_Circ VS My program
    //                 CIRI
    //                 Circ_Explorer
    cout << endl << endl;
    cout << "************ Find_Circ Compare with My Result, CIRI, and Circ_Explorer ************" << endl;
    cout << "************ ********************************************************* ************" << endl;
    pResultComparison->CheckIntersetBTMyProgramAndFindCirc(stConfig.strCurHit, stConfig.strFindCircHit);
    pResultComparison->CheckIntersetBTCIRIAndFindCirc(stConfig.strCIRIHit, stConfig.strFindCircHit);
    pResultComparison->CheckIntersetBTCircExplorerAndFindCirc(stConfig.strCircExplorerHit, stConfig.strFindCircHit);

    // 比对 Circ_Explorer VS My program
    //                      CIRI
    //                      Find_Circ
    cout << endl << endl;
    cout << "************ Circ_Explorer Compare with My Result, CIRI, and Find_Circ ************" << endl;
    cout << "************ ********************************************************* ************" << endl;
    pResultComparison->CheckIntersetBTMyProgramAndCircExplorer(stConfig.strCurHit, stConfig.strCircExplorerHit);
    pResultComparison->CheckIntersetBTCIRIAndCircExplorer(stConfig.strCIRIHit, stConfig.strCircExplorerHit);
    pResultComparison->CheckIntersetBTFindCircAndCircExplorer(stConfig.strFindCircHit, stConfig.strCircExplorerHit);

#endif

    delete pResultComparison;
    pResultComparison = NULL;

    ClsTroubleShoot* pTS = new ClsTroubleShoot();
    //pTS->CheckBWAMappingResult(stConfig.strBWA);
    delete pTS;
    pTS = NULL;
//#endif

    //Release resources
    delete pKT;
    pKT = NULL;
    delete pFindCandi;
    pFindCandi = NULL;

    return 0;
}

struct St_FindCircKT
{
    St_Fasta* pRef;
    ClsKmerTable* pKT;
    vector<St_Fastq>* pReads;
    St_Row_Chrom* pRowChrom;
    vector<St_Row_Chrom>* pWholeChrom;
    ClsFindCandidate* pFindCandi;
    St_Config* pConfig;
    string strChromName;
    int iChromIndex;

    St_FindCircKT()
    {
        Reset();
    }

    void Init(St_Fasta* pV1, ClsKmerTable* pV2, vector<St_Fastq>* pV3,
              St_Row_Chrom* pV4, vector<St_Row_Chrom>* pV5,
              ClsFindCandidate* pV6, St_Config* pV7, string strV8, int iV9)
    {
        pRef = pV1;
        pKT = pV2;
        pReads = pV3;
        pRowChrom = pV4;
        pWholeChrom = pV5;
        pFindCandi = pV6;
        pConfig = pV7;
        strChromName = strV8;
        iChromIndex = iV9;
    }

    void Reset()
    {
        pRef = NULL;
        pKT = NULL;
        pReads = NULL;
        pRowChrom = NULL;
        pWholeChrom = NULL;
        pFindCandi = NULL;
        pConfig = NULL;
        strChromName = "";
        iChromIndex = -1;
    }
};

void* FindCircForSingleChrom(void* pValue)
{
    cout << "main() : creating thread " << endl;

    St_FindCircKT* pFindCircKT = (St_FindCircKT*)pValue;

    map<unsigned int, vector<St_PosInfo> > mpKT;
    pFindCircKT->pKT->CreateKmerTable( mpKT,
                                       pFindCircKT->pConfig->strRefPath, pFindCircKT->pConfig->iReadsLen,
                                       pFindCircKT->pConfig->fKmerRatio,
                                       pFindCircKT->pRowChrom, pFindCircKT->pRef, pFindCircKT->iChromIndex);

    pFindCircKT->pFindCandi->CheckHitting(pFindCircKT->pConfig->iMinSupportReads,
                                          pFindCircKT->pConfig->fKmerRatio, pFindCircKT->pConfig->iReadsLen,
                                          mpKT, *(pFindCircKT->pWholeChrom),
                                          *(pFindCircKT->pReads), pFindCircKT->strChromName);

    pthread_exit(NULL);
}

void FindCircRNA(St_Config& stConfig, ClsKmerTable* pKT, ClsFindCandidate* pFindCandi,
                 vector<St_Fastq>& vFastq, vector<St_Row_Chrom>& vChrom, vector<St_Fasta>& vFasta)
{

    cout << "The second vChrom Size: " << vChrom.size() << endl;

    cout << "Step 3: Detect Circular RNA" << endl;
    string strCmd = (string)"mkdir -p " + get_current_dir_name() + "/Detection_Result";
    system(strCmd.c_str());
    strCmd = "rm ./Detection_Result/*";
    system(strCmd.c_str());

    int iThreadNum = 0;
    vector<St_FindCircKT> vFindCircKT;
    St_FindCircKT stFindCircKT;
    for(vector<St_Fasta>::iterator itrRef = vFasta.begin(); itrRef != vFasta.end(); itrRef++)
    {
        //cout << itrRef->strName << endl;
        string strRefName = "";
        if(itrRef->strName.find(' ') == string::npos)
            strRefName = itrRef->strName;
        else
            strRefName = itrRef->strName.substr(0, itrRef->strName.find(' '));

        int iChromIndex = 0;
        for(vector<St_Row_Chrom>::iterator itrChrom = vChrom.begin(); itrChrom != vChrom.end(); itrChrom++)
        {
            if(strRefName == itrChrom->strName)
            {
//                map<unsigned int, vector<St_PosInfo> > mpKT;
//                pKT->CreateKmerTable(mpKT, stConfig.strRefPath, stConfig.iReadsLen, stConfig.fKmerRatio,
//                                     &(*itrChrom), &(*itrRef), iChromIndex);

//                //1: Collect Hit Time and Info Order
//                cout << "Collect Hit Time and Info Order" << endl;
//                cout << "mpKT size: " << mpKT.size() << endl;
//                cout << "We got size!" << endl;

//                int i = 0;
//                for(map<unsigned int, vector<St_PosInfo> >::iterator itr = mpKT.begin(); itr != mpKT.end(); itr++)
//                {
//                    if(i > 10)
//                        break;

//                    cout << "---" << endl;
//                    cout << itr->first << endl;
//                    cout << IntToStr(itr->second.begin()->ucChromIndex)
//                         << " --- " << IntToStr(itr->second.begin()->ucTranscriptIndex) << endl;
//                    cout << "---" << endl;

//                    i++;
//                }

                //-->Init St_FindCircKT
                stFindCircKT.Reset();
                stFindCircKT.Init(&(*itrRef), pKT, &vFastq, &(*itrChrom), &vChrom,
                                  pFindCandi, &stConfig, strRefName, iChromIndex);
                vFindCircKT.push_back(stFindCircKT);
                iThreadNum++;
                break;
            }
            else
                iChromIndex++;
        }
    }

    cout << "The third vChrom Size: " << vChrom.size() << endl;

    //Step 1: Define how may threads you want to use
    pthread_t threads[iThreadNum];
    pthread_attr_t attr;

    //Step 2: Initialize and set thread joinable
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    //Step 3: Do the main body of multiple threads function
    cout << "Thread Num: " << iThreadNum << endl;
    for(int i = 0; i < iThreadNum; i++)
    {
        //cout << "main() : creating thread, " << i << endl;

        int rc = pthread_create(&threads[i], NULL, FindCircForSingleChrom, (void *)&vFindCircKT[i]);

        if(rc)
        {
            cout << "Error:unable to create thread," << rc << endl;
            exit(-1);
        }
    }

    //4: free attribute and wait for the other threads
    pthread_attr_destroy(&attr);

    void *status;
    //5: Join those threads together --> to make sure the remaining part of the code will be run only after those threads been finished
    for(int i = 0; i < iThreadNum; i++)
    {
        int rc = pthread_join(threads[i], &status);
        if (rc)
        {
            cout << "Error:unable to join," << rc << endl;
            exit(-1);
        }

        //cout << "Main: completed thread id :" << i ;
        //cout << "  exiting with status :" << status << endl;
    }

    //6: Do the remaning thing
    cout << endl << "ALL SET!!!" << endl;

    cout << "Combine the result together" << endl;
    strCmd = "rm ./Detection_Result/Brief_sum.txt";
    system(strCmd.c_str());
    strCmd = "cat ./Detection_Result/Brief* > ./Detection_Result/Brief_sum.txt";
    system(strCmd.c_str());
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Test(string strFaPath)
{

    St_Candidate stCurCandi1(0, 1, 1);
    St_Candidate stCurCandi2(0, 2, 2);
    map<St_Candidate, int> mpvCandi;
    mpvCandi[stCurCandi1] = 1;
    mpvCandi[stCurCandi2] = 1;
    cout << IntToStr(mpvCandi[stCurCandi1]) << endl;
    cout << IntToStr(mpvCandi[stCurCandi2]) << endl;
    if(mpvCandi.find(stCurCandi2) != mpvCandi.end())
        mpvCandi[stCurCandi2]++;
    cout << IntToStr(mpvCandi[stCurCandi2]) << endl;


    cout << GetReverseCompelement("CCCAATGGCTGAGGGTGCAGACCATTCTCAATTCCTATTTCAAATCCTCCACCAAGGCCAGCTATCCCAATGGCTGAGGGTGCAGACCATTCTCAATTCCTATTTCAAATCCTCCACCAAGGCCAGCTATCCCAATGGCTGAGGGTGCAG")
         << endl;


    ClsFastaReader* pFaReader = new ClsFastaReader();
    vector<St_Fasta> vFa;
    pFaReader->ReadFastaRegular(strFaPath, vFa);

    cout << vFa[2].strSeq.substr(383816, 71) << endl << endl;
    //cout << vFa[2].strSeq.substr(382474, 198) << endl;

    delete pFaReader;
    pFaReader = NULL;

    map<St_PosInfo, int> mpInfo;

    St_PosInfo stInfo2(1,1,1,1);
    St_PosInfo stInfo1(2,2,2,2);
    mpInfo[stInfo1] = 2;
    mpInfo[stInfo2] = 1;

    string strTmp = "AAAAAA";
    if(strTmp.find('N') != string::npos)
    {
        cout << "Got it" << endl;
    }

    if(strTmp.find('N') == string::npos)
    {
        cout << "Cannot fint N" << endl;
    }

    unsigned char cT = 2;

    if(cT == 2)
        cout << "I am 2" << endl;
    else
        cout << "Hi" << endl;
    //int i = 0;
}

