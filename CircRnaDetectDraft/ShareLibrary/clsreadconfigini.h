#ifndef CLSREADCONFIGINI_H
#define CLSREADCONFIGINI_H
#include "string"
#include <map>
#include <vector>
using namespace std;

/*********************************************
 * How to use it (example):
 *
 * 1: Identify the structure of your config file
 *  struct St_Config
    {
        string strRef;
        string strOrgLR; //LR means long reads
        string strCorrectedLR;
        string strOrgSR; //SR means short reads
    };

   2: Read the the value of each variants in your config file
      Notice: the first parameter is the file path of your config file

    St_Config stConfig;
    void ReadConfig(St_Config& stConfig, int argc, char **argv)
    {
        //Reads the disk file
        St_Config stConfig;
        ClsReadConfigIni* pIni = new ClsReadConfigIni();
        pIni->ReadIni(argv[1]); // The first valid parameters
        for(vector<St_Section>::iterator itr = pIni->GetConfigInfo().vSect.begin();
            itr != pIni->GetConfigInfo().vSect.end(); itr++)
        {
            if(itr->strName == "General")
            {
                for(map<string, string>::iterator itrmp = itr->m_mpKeyValue.begin();
                    itrmp != itr->m_mpKeyValue.end(); itrmp++)
                {
                    if(itrmp->first == "Reference")
                        stConfig.strRef = itrmp->second;
                    else if(itrmp->first == "OriginalLongReads")
                        stConfig.strOrgLR = itrmp->second;
                    else if(itrmp->first == "CorrectedLongReads")
                        stConfig.strCorrectedLR = itrmp->second;
                    else if(itrmp->first == "OriginalShortReads")
                        stConfig.strOrgSR = itrmp->second;
                }
            }
        }
        delete pIni;
        pIni = NULL;
    }
 *******************************************************************
 */

struct St_Section
{
    string strName;
    map<string, string> m_mpKeyValue;

    St_Section():strName("")
    {}

    void Init()
    {
        strName = "";
        m_mpKeyValue.clear();
    }
};

struct St_ConfigIni
{
    vector<St_Section> vSect;
};

class ClsReadConfigIni
{
public:
    ClsReadConfigIni();
    ~ClsReadConfigIni();

public:
    void ReadIni(char* czIniPath);

    St_ConfigIni& GetConfigInfo()
    {
        return m_stIni;
    }

private:
    St_ConfigIni m_stIni;
};

#endif // CLSREADCONFIGINI_H
