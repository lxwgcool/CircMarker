#include <iostream>
#include "../../ShareLibrary/clsgtfparse.h"
#include "clsconfig.h"
#include "../../ShareLibrary/clskmeralgorithm.h"
#include "../../ShareLibrary/clsbasealgorithm.h"
#include "clskmertable.h"
#include "clsfindcandidate.h"
#include "clsresultcomparison.h"
#include "clstroubleshoot.h"

using namespace std;

void Test(string strFaPath);

#define ONLY_DO_COMPARISON
//#define CUR_COMPARISON
//#define SIMULATION_COMPARISON
#define CHECK_INTERSECTION

int main(int argc, char **argv)
{   
    //1: Read Ini
    ClsConfig* pConfig = new ClsConfig();
    St_Config stConfig;
    pConfig->ReadConfig(stConfig, argv[1]);
    delete pConfig;
    pConfig = NULL;

    //Test(stConfig.strRefPath);

    //2: Read GTF---------------------------
    //   Get the information: Chromosn, Gene, Transcript, and Exon
    vector<St_Row_Chrom> vChrom;
    ClsGTFParse* pGTFParse = new ClsGTFParse();
    pGTFParse->ReadGTF(stConfig.strGtfPath, vChrom); //Done
//#ifndef ONLY_DO_COMPARISON
    pGTFParse->GetTagValue(stConfig.strRefPath, stConfig.iReadsLen, stConfig.fKmerRatio, KMERLEN,
                           vChrom);
//#endif
    delete pGTFParse;
    pGTFParse = NULL;

    //3: Create Kmer Table-----------------
    ClsKmerTable* pKT = new ClsKmerTable();
#ifndef ONLY_DO_COMPARISON
    cout << endl << "--------------CreateKmerTable--------------" << endl << endl;
    pKT->CreateKmerTable(stConfig.strRefPath, stConfig.iReadsLen, stConfig.fKmerRatio, vChrom);
#endif

    //4: Check Hitting
    ClsFindCandidate* pFindCandi = new ClsFindCandidate();
#ifndef ONLY_DO_COMPARISON
    pFindCandi->CheckHitting(stConfig.strReads1Path, stConfig.strReads2Path,
                             stConfig.iMinSupportReads, stConfig.fKmerRatio,
                             pKT->GetKT(), vChrom);
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
        cout << "Fuck" << endl;
    int i = 0;
}

