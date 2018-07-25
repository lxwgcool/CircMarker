#ifndef CLSRESULTCOMPARISON_H
#define CLSRESULTCOMPARISON_H
#include "clsfindcandidate.h"

class ClsResultComparison
{
public:
    ClsResultComparison();

public:
    //Compare with simulation: --->
    void CompareMyResultWithSimulation(string strCurPath, string strSimulationResult, int iMaxSupportReads);
    void CompareCiriWithSimulation(string strCIRIResult, string strSimulationResult);
    void CompareFindCircWithSimulation(string strFindCircResult, string strSimulationResult);
    void CompareCircExplorerWithSimulation(string strCIRCexplorerResult, string strSimulationResult);
    //<---

    //Check intersection --->
    void CheckIntersectionOfMyProgram(string strCurResult,
                                      string strRibominusCur, int iMaxSupportReads); // Between ribominus data and RNaseR & PolyA- data
                                                                              // LC means: linear and Circular Path
                                                                              // COnly means: Circular Only
    void CheckIntersectionOfCIRI(string strCIRIResult, string strRibominusCIRI); //  ...
    void CheckIntersectionOfFindCirc(string strFindCircResult, string strRibominusFindCirc); //  ...
    void CheckIntersectionOfCircExplorer(string strCIRCexplorerResult, string strRibominusCIRCexplorer); //  ...
    void GetMostSupportedCandiAndCheckIntersection(string strCurResult,
                                                   string strCIRIResult,
                                                   string strFindCircResult,
                                                   string strCIRCexplorerResult);
    void EraseDuplicate(vector<St_Candidate>& vResult);
    //<---

    void CompareGTFWithCircBaseResult(vector<St_Row_Chrom>& vChrom, string strCircBasePath);

    //Compare result with CircBaseResult, include: My Result, CIRI, Find_Circ, CircExplorer
    void CompareMyResultWithCircBaseResult(vector<St_Row_Chrom>& vChrom, string strCircBasePath,
                                           string strCurResultPath, int iCurChromIndex=-1);

    void CompareCiriWithCircBaseResult(vector<St_Row_Chrom>& vChrom, string strCircBasePath,
                                       string strCiriResultPath, int iCurChromIndex=-1);

    void CompareFindCircWithCircBaseResult( vector<St_Row_Chrom>& vChrom, string strCircBasePath,
                                            string strFindCircPath, int iCurChromIndex=-1);

    void CompareCircExplorerWithCircBaseResult( vector<St_Row_Chrom>& vChrom,
                                                string strCircBasePath,
                                                string strCircExplorerPath, int iCurChromIndex=-1);
    void SummaryCompareResult();

    // 比对 My program VS CIRI
    //                   Find_Circ
    //                   Circ_Explorer
    void CheckIntersetBTCIRIAndMyProgram(string strCIRI, string strMyProgram);
    void CheckIntersetBTFindCircAndMyProgram(string strFindCirc, string strMyProgram);
    void CheckIntersetBTCircExplorerAndMyProgram(string strCircExplorer, string strMyProgram);

    // 比对 CIRI VS My program
    //             Find_Circ
    //             Circ_Explorer
    void CheckIntersetBTMyProgramAndCIRI(string strMyProgram, string strCIRI);
    void CheckIntersetBTFindCircAndCIRI(string strFindCirc, string strCIRI);
    void CheckIntersetBTCircExplorerAndCIRI(string strCircExplorer, string strCIRI);

    // 比对 Find_Circ VS My program
    //                 CIRI
    //                 Circ_Explorer
    void CheckIntersetBTMyProgramAndFindCirc(string strMyProgram, string strFindCirc);
    void CheckIntersetBTCIRIAndFindCirc(string strCIRI, string strFindCirc);
    void CheckIntersetBTCircExplorerAndFindCirc(string strCircExplorer, string strFindCirc);

    // 比对 Circ_Explorer VS My program
    //                      CIRI
    //                      Find_Circ
    void CheckIntersetBTMyProgramAndCircExplorer(string strMyProgram, string strCircExplorer);
    void CheckIntersetBTCIRIAndCircExplorer(string strFindCirc, string strCircExplorer);
    void CheckIntersetBTFindCircAndCircExplorer(string strCIRI, string strCircExplorer);

private:
    // Brief.txt (same directory as executable file)
    void ParseBrief(string strCurPath, vector<St_Candidate>& vCandi); // Brief is my current result
                                                                      // It is also the common format of comparison
    // Simulation Result
    void ParseSimulationResult(string strStdPath, vector<St_Candidate>& vCandi);

    // CIRI Result
    void ParseCIRIResult(string strCIRIPath, vector<St_Candidate>& vCandi);

    //CIRCexplorer Result
    void ParseCircExplorerResult(string strCIRCexplorerPath, vector<St_Candidate>& vCandi);

    //Find_Circ Result
    void ParseFindCircResult(string strFindCircPath, vector<St_Candidate>& vCandi);

    //Parse CircBase Std Circular RNA
    void ParseCircBaseResult(string strCircBasePath, vector<St_Candidate>& vCandi);

private: //sub function 比如一些类似的获取两边是否匹配的结论等
    int CheckBothHitExonBoundary(vector<St_Row_Chrom>& vChrom, St_Candidate* pCandi);
    bool CompareTwoCandi(St_Candidate* pCandi1, St_Candidate* pCandi2);

private:
    //Record some variables
    int iMyProgramHitNum;
    int iCIRIHitNum;
    int iFindCircHitNUm;
    int iCircExplorerHitNum;

    float fMyProgramTP;
    float fCIRITP;
    float fFindCircTP;
    float fCircExplorerTP;
};

#endif // CLSRESULTCOMPARISON_H
