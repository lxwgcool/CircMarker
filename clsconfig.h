#ifndef CLSCONFIG_H
#define CLSCONFIG_H

#include <string>
using namespace std;

struct St_Config
{
    //General
    string strRefPath;
    string strGtfPath;
    string strReads1Path;
    string strReads2Path;

    //Parameters
    int    iMinSupportReads;
    int    iMaxSupportReads;
    int    iReadsLen;
    float  fKmerRatio; // what's the length ratio of reads length
                       //(总kmer的数量for one exon is 2*(iKmerRatio% * Len(reads))

    //Result
    string strCurResult;
    string strSimulationResult;
    string strCIRIResult;
    string strCIRCexplorerResult;
    string strFindCircResult;
    string strCircBaseResult;    

    string strRibominusCur;
    string strRibominusCIRI;
    string strRibominusFindCirc;
    string strRibominusCIRCexplorer;

    //Hitting Comparison
    string strCIRIHit;
    string strCurHit;
    string strFindCircHit;
    string strCircExplorerHit;
    int iCurChromIndex;

    //Mapping
    string strBWA;

    St_Config():strRefPath(""), strGtfPath(""), strReads1Path(""), strReads2Path(""),
                iMinSupportReads(1), iMaxSupportReads(999), iReadsLen(101), fKmerRatio(.3),
                strCurResult(""), strSimulationResult(""), strCIRIResult(""), strCIRCexplorerResult(""),
                strFindCircResult(""), strCircBaseResult(""), strRibominusCur(""), strRibominusCIRI(""),
                strRibominusFindCirc(""), strRibominusCIRCexplorer(""),strCIRIHit(""), strCurHit(""),
                strFindCircHit(""), strCircExplorerHit(""), iCurChromIndex(-1), strBWA("")
    {}
};

class ClsConfig
{
public:
    ClsConfig();

public:
    void ReadConfig(St_Config& stConfig, char* cpIniPath);
};

#endif // CLSCONFIG_H
