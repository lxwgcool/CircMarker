#include "clsreadconfigini.h"
#include "unistd.h"
#include <fstream>
#include <iostream>

ClsReadConfigIni::ClsReadConfigIni()
{
}

ClsReadConfigIni::~ClsReadConfigIni()
{}

void ClsReadConfigIni::ReadIni(char* czIniPath)
{
    //Check if such file existed
    if(::access(czIniPath, 0) != 0)
    {
        cout << "Ini File Do not Existed!" << endl;
        return;
    }
    m_stIni.vSect.clear();
    St_Section stSection;
    ifstream ifs;
    ifs.open(czIniPath);
    string strLine = "";
    while(!ifs.eof())
    {
        //Get the value of  current line
        getline(ifs, strLine);
        //For the last abnormal item
        if(ifs.eof())
        {
            if(stSection.strName != "")
                m_stIni.vSect.push_back(stSection);
            break;
        }
        //For the normal item
        if(strLine == "" ||  // Null
           strLine == "\n" || // Nil Enter
           strLine.at(0) == '#') // Comments
            continue;
        //For the start of section
        if(strLine.find('[') != string::npos && strLine.find(']') != string::npos) // The new section
        {
            //Record the old section into
            if(stSection.strName != "")
                m_stIni.vSect.push_back(stSection);
            //Init the new section
            stSection.Init();
            string::size_type iStart = strLine.find('[') + 1;
            string::size_type iEnd = strLine.find(']');
            stSection.strName = strLine.substr(iStart, iEnd - iStart);
        }
        else //2: Extract the value item of current line
        {
            string::size_type iSplit = strLine.find('=');
            if(iSplit == string::npos)
                continue;
            stSection.m_mpKeyValue[strLine.substr(0, iSplit)] = strLine.substr(iSplit+1,
                                                                               strLine.length() - iSplit - 1);
        }
         // for the last item
        if(ifs.eof())
        {
            if(stSection.strName != "")
                m_stIni.vSect.push_back(stSection);
        }
    }
}

