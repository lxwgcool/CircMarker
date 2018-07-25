#include "clsconfig.h"
#include "../../ShareLibrary/clsreadconfigini.h"
#include "stdlib.h"
#include <iostream>
using namespace std;

ClsConfig::ClsConfig()
{
}

void ClsConfig::ReadConfig(St_Config& stConfig, char* cpIniPath)
{
    ClsReadConfigIni* pIni = new ClsReadConfigIni();
    pIni->ReadIni(cpIniPath); // The first valid parameters

    for(vector<St_Section>::iterator itr = pIni->GetConfigInfo().vSect.begin();
        itr != pIni->GetConfigInfo().vSect.end(); itr++)
    {
        if(itr->strName == "General")
        {
            for(map<string, string>::iterator itrmp = itr->m_mpKeyValue.begin();
                itrmp != itr->m_mpKeyValue.end(); itrmp++)
            {
                if(itrmp->first == "Reference")
                    stConfig.strRefPath = itrmp->second;
                else if(itrmp->first == "GTF")
                    stConfig.strGtfPath = itrmp->second;
                else if(itrmp->first == "Reads1")
                    stConfig.strReads1Path = itrmp->second;
                else if(itrmp->first == "Reads2")
                    stConfig.strReads2Path = itrmp->second;                
            }
        }
        else if(itr->strName == "Parameters")
        {
            for(map<string, string>::iterator itrmp = itr->m_mpKeyValue.begin();
                itrmp != itr->m_mpKeyValue.end(); itrmp++)
            {
                if(itrmp->first == "MinSupportReads")
                    stConfig.iMinSupportReads = atoi(itrmp->second.c_str());
                if(itrmp->first == "MaxSupportReads")
                    stConfig.iMaxSupportReads = atoi(itrmp->second.c_str());
                else if(itrmp->first == "ReadsLen")
                    stConfig.iReadsLen = atoi(itrmp->second.c_str());
                else if(itrmp->first == "KmerRatio")
                    stConfig.fKmerRatio = (float)atoi(itrmp->second.c_str()) / 100;
            }
        }
        else if(itr->strName == "Results")
        {
            for(map<string, string>::iterator itrmp = itr->m_mpKeyValue.begin();
                itrmp != itr->m_mpKeyValue.end(); itrmp++)
            {
                if(itrmp->first == "CurResult")
                    stConfig.strCurResult = itrmp->second;
                else if(itrmp->first == "SimulationResult")
                    stConfig.strSimulationResult = itrmp->second;
                else if(itrmp->first == "CIRIResult")
                    stConfig.strCIRIResult = itrmp->second;
                else if(itrmp->first == "CIRCexplorerResult")
                    stConfig.strCIRCexplorerResult = itrmp->second;
                else if(itrmp->first == "FindCircResult")
                    stConfig.strFindCircResult = itrmp->second;
                else if(itrmp->first == "CircBaseResult")
                    stConfig.strCircBaseResult = itrmp->second;
                else if(itrmp->first == "RibominusCur")  //--------> For ribominus comparison
                    stConfig.strRibominusCur = itrmp->second;
                else if(itrmp->first == "RibominusCIRI")
                    stConfig.strRibominusCIRI = itrmp->second;
                else if(itrmp->first == "RibominusFindCirc")
                    stConfig.strRibominusFindCirc = itrmp->second;
                else if(itrmp->first == "RibominusCIRCexplorer")
                    stConfig.strRibominusCIRCexplorer = itrmp->second; //<-----
            }
        }
        else if(itr->strName == "HitCompare")
        {
            for(map<string, string>::iterator itrmp = itr->m_mpKeyValue.begin();
                itrmp != itr->m_mpKeyValue.end(); itrmp++)
            {
                if(itrmp->first == "CIRIHit")
                    stConfig.strCIRIHit = itrmp->second;
                else if(itrmp->first == "CurHit")
                    stConfig.strCurHit = itrmp->second;
                else if(itrmp->first == "FindCircHit")
                    stConfig.strFindCircHit = itrmp->second;
                else if(itrmp->first == "CircExplorerHit")
                    stConfig.strCircExplorerHit = itrmp->second;
                else if(itrmp->first == "CurChromIndex")
                    stConfig.iCurChromIndex = atoi(itrmp->second.c_str());
            }
        }
        else if(itr->strName == "Mapping")
        {
            for(map<string, string>::iterator itrmp = itr->m_mpKeyValue.begin();
                itrmp != itr->m_mpKeyValue.end(); itrmp++)
            {
                if(itrmp->first == "BWA")
                    stConfig.strBWA = itrmp->second;
            }
        }
    }

    delete pIni;
    pIni = NULL;
}

bool ClsConfig::CheckConfig(St_Config& stConfig)
{
    bool bReadsPath = true;
    bool bReferencePath = true;
    bool bAnnotationPath = true;
    if(stConfig.strReads1Path == "") // do not care about reads2path
    {
        cout << "Sequence Reads is mandated" << endl;
        bReadsPath = false;
    }

    if(stConfig.strRefPath == "")
    {
        cout << "Reference file is mandated!" << endl;
        bReferencePath = false;
    }

    if(stConfig.strGtfPath == "")
    {
        cout << "Annotation File is mandated!" << endl;
        bAnnotationPath = false;
    }

    if(!bReadsPath || !bReferencePath || !bAnnotationPath)
    {
        cout << "**********************************" << endl;
        cout << "How to use CircMarker" << endl;
        cout << "./CircRnaDetectDraft ./config.ini" << endl;
        cout << "**********************************" << endl;
        return false;
    }
    else
        return true;
}


