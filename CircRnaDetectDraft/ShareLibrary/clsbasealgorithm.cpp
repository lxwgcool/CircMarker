/******************This is cpp file (*.cpp) of basical function**************
 * *************************************************************************/
#include "clsbasealgorithm.h"
#include <memory>
#include <functional>
#include "string.h" // For memset
#include <stdio.h>
#include <cctype>
#include <iostream>
#include <algorithm>
#include <vector>
#include "unistd.h" // For "getcwd"
#include <math.h>
#include <fstream>
#include <dirent.h> // For the stucture "DIR"
#include <sys/stat.h>
#define MAXBUFSIZE 260

string& ToUpper(string& strValue)
{
    if(strValue == "")
        return strValue;
    transform(strValue.begin(), strValue.end(), strValue.begin(), (int(*)(int))toupper);
    return strValue;
}

string& ToLower(string& strValue)
{
    if(strValue == "")
        return strValue;
    transform(strValue.begin(), strValue.end(), strValue.begin(), (int(*)(int))tolower);
    return strValue;
}

string IntToStr(int iValue)
{
    char czStr[20];
    ::memset(czStr, 0, 20);
    ::sprintf(czStr, "%d", iValue);
    return czStr;
}

string FloatToStr(double fValue, int iDigit)
{
    char czStr[20];
    ::memset(czStr, 0, 20);

    ::sprintf(czStr, ("%." + IntToStr(iDigit) + "f").c_str(), fValue);
    return czStr;
}

using namespace std;

string& ltrim(string &str) {
    string::iterator p = find_if(str.begin(), str.end(), not1(ptr_fun<int, int>(isspace)));
    str.erase(str.begin(), p);
    return str;
}

string& rtrim(string &str) {
    string::reverse_iterator p = find_if(str.rbegin(), str.rend(), not1(ptr_fun<int , int>(isspace)));
    str.erase(p.base(), str.end());
    return str;
}

string& trim(string &str) {
    ltrim(rtrim(str));
    return str;
}

string GetReverse(string strOrg)
{
    //wanna the reverse complementary sequence
    //1: Get reverse
    reverse(strOrg.begin(), strOrg.end());
    return strOrg;
}

string GetHigherFolderPath(string strCurPath, int iLevel)
{
    vector<int> vLevel;
    size_t iPos = 0;
    while((iPos = strCurPath.find('/', iPos)) != string::npos)
    {
        vLevel.push_back(iPos);
        iPos += 1;
    }
    if((int)vLevel.size() < iLevel)
        return "";
    if(iLevel == 0)
        return  strCurPath.at(strCurPath.length()-1) != '/' ? strCurPath + "/" : strCurPath;
    return strCurPath.substr(0, vLevel[vLevel.size()-iLevel]+1);
}

string GetCurExeFolderPath()
{
    //::get_current_dir_name();
    char buf[MAXBUFSIZE];
    memset(buf, 0, MAXBUFSIZE);
    ::getcwd(buf, MAXBUFSIZE);
    return string(buf);
}

bool IsMissing(char nt)
{
    return nt == 'N' || nt == 'n';
}

char GetComplement(char bp)
{
    //
    if( IsMissing(bp) )
    {
        return bp;
    }
    //char bpUse = toupper(bp);
    if( bp == 'A')
    {
        return 'T';
    }
    else if(bp == 'a')
    {
        return 't';
    }
    else if( bp == 'T')
    {
        return 'A';
    }
    else if( bp == 't')
    {
        return 'a';
    }
    else if( bp == 'G')
    {
        return 'C';
    }
    else if(bp == 'g')
    {
        return 'c';
    }
    else if(bp == 'C')
    {
        return 'G';
    }
    else if(bp == 'c')
    {
        return 'g';
    }
    return 'N';
}

string GetReverseCompelement(string strOrg, bool bRevsCompelement)
{
    if(!bRevsCompelement) //Do not need the reverse complementary sequence
        return strOrg;
    //wanna the reverse complementary sequence
    //1: Get reverse
    reverse(strOrg.begin(), strOrg.end());
    //2: Get Compement
    for(unsigned int i=0; i<strOrg.length(); i++)
    {
        strOrg[i] = GetComplement(strOrg[i]);
    }
    return strOrg;
}

//Display the string in file --> multi-lines
void DisplayString(ofstream& ofs, string& strValue, int iLenPerLine)
{
    if(strValue == "")
        return;

    int iLineNum = ceil((float)strValue.length() / iLenPerLine);
    for(int i=0; i<iLineNum; i++)
    {
        if(i + 1 == iLineNum) // This is the last item
            ofs << strValue.substr(i*iLenPerLine, strValue.length() - i*iLenPerLine) << endl;
        else //This is not the last item
            ofs << strValue.substr(i*iLenPerLine, iLenPerLine) << endl;
    }
}

//Get N,n pos
string::size_type GetNextNPos(string& strOrg, int iCurPos, bool bNormalOrientation)
{
    string::size_type sztPos = string::npos;
    if(bNormalOrientation) //正常的秩序
    {
        sztPos = strOrg.find("N", iCurPos);
        if(sztPos == string::npos)
            sztPos = strOrg.find("n", iCurPos);
    }
    else //The reverse orientation
    {
        sztPos = strOrg.rfind("N", iCurPos);
        if(sztPos == string::npos)
            sztPos = strOrg.rfind("n", iCurPos);
    }
    return sztPos;
}

bool CheckNPos(char cCurChar)
{
    if( cCurChar == 'N' || cCurChar == 'n')
        return true;
    else
        return false;
}

string GetStringItemValue(string& strValue, string::size_type sztpPos, bool bForward, char cSplit)
{
    string strItem = ""; // notice do NOT include the character in the position of "sztpPos"
    if(bForward) //Forward direction 正向序列
    {
        string::size_type sztpHit = strValue.find(cSplit, sztpPos + 1);
        if(sztpHit == string::npos)
            return "";
        else
        {
            int iLen = sztpHit - sztpPos - 1;
            if(iLen > 0)
            {
                strItem = strValue.substr(sztpPos + 1, iLen);
                return strItem;
            }
            else
                return "";
        }
    }
    else //Backward direction 反向序列
    {
        string::size_type sztpHit = strValue.rfind(cSplit, sztpPos - 1);
        if(sztpHit == string::npos)
            return "";
        else
        {
            int iLen = sztpPos - sztpHit - 1;
            if(iLen > 0)
            {
                strItem = strValue.substr(sztpHit + 1, iLen);
                return strItem;
            }
            else
                return "";
        }
    }
}

string GetRatio(float fValue, int iPrecision, bool bPrintSymbol)
{
    char czValue[16];
    memset(czValue, 0, 16);    
    if(bPrintSymbol)
    {
        string strFormate = "%." + IntToStr(iPrecision) + "f";
        sprintf(czValue, strFormate.c_str(), fValue*100);
        return czValue + (string)"%";
    }
    else // the dot is based on value%.
    {
        string strFormate = "%." + IntToStr(iPrecision + 2) + "f";
        sprintf(czValue, strFormate.c_str(), fValue);
        return czValue;
    }
}

//czFolderPath --> 搜索的目录
void GetAllFilesInFolder(vector<string>& vFileName, const char* czFolderPath)
{
    //clear old container
    vFileName.clear();

    //parse the file path in folder
    int return_code;
    DIR *dir;
    struct dirent entry;
    struct dirent *res;
    if ((dir = opendir(czFolderPath)) != NULL)
    {//打开目录
       for(return_code = readdir_r(dir, &entry, &res);
           res != NULL && return_code == 0;
           return_code = readdir_r(dir, &entry, &res))
       {
           if(entry.d_type != DT_DIR)
           {//存放到列表中               
               vFileName.push_back(string(entry.d_name));
           }
       }
       closedir(dir);//关闭目录
    }
}

string GetCompleCode(int iOrg, int iStdDigit)
{
    string strCode = "";
    int iOrgDigitNum = 0;
    if(iOrg == 0)
        iOrgDigitNum = 1;
    else
    {
        int iTempOrg = iOrg;
        while(iTempOrg)
        {
            iTempOrg /= 10;
            iOrgDigitNum++;
        }
    }
    if(iOrgDigitNum >= iStdDigit)
        strCode = ::IntToStr(iOrg);
    else
    {
        for(int i=0; i< iStdDigit - iOrgDigitNum; i++)
        {
            strCode += "0";
        }
        strCode += ::IntToStr(iOrg);
    }
    return strCode;
}

string GetHMSTimeFormat(double dSumSeconds)
{
    int iHour, iMinute, iSecond;
    iHour = dSumSeconds / (60 * 60);
    iMinute = (dSumSeconds - (iHour * 60 * 60)) / 60;
    iSecond = dSumSeconds - (iHour * 60 * 60) - (iMinute * 60);
    string strHMSFormat = IntToStr(iHour) + "h " + IntToStr(iMinute) + "m " +
                          IntToStr(iSecond) + "s";
    return strHMSFormat;
}

bool CheckContainState(int iStart, int iEnd, int iPos,
                       bool bTolerant, int iOffSet)
{
    bool bContain = false;
    if(bTolerant)
    {
        if(iPos > (iStart - iOffSet) && iPos < (iEnd + iOffSet)) // offset是对查询范围的左右offset距离的扩展
            bContain = true;
        else if(iPos < (iStart + iOffSet) && iPos > (iEnd - iOffSet)) // 这意味着Start和End位置互调了
            bContain = true;
    }
    else
    {
        if(iPos > iStart && iPos < iEnd)
            bContain = true;
        else if(iPos < iStart && iPos > iEnd) // 这意味着Start和End位置互调了
            bContain = true;
    }
    return bContain;
}

/*
int GetMatchPartLen(BamAlignment& al)
{
    int iMatchLen = 0;
    for(std::vector<CigarOp>::iterator itr = al.CigarData.begin();
        itr != al.CigarData.end(); itr++)
    {
        switch(itr->Type)
        {
            case 'M': // alignment match (can be a sequence match or mismatch)
                iMatchLen += itr->Length;
                break;
            default:
                break;
        }
    }
    return iMatchLen;
}*/

/*
string GetSoftClipPart(BamAlignment& al, En_ClipPart enPart) // The part that been cut off
{    
    //目前的算法是可以求得第一次出现softclip
    int iStart = 0;
    vector<string> vClispSeq;
    for(std::vector<CigarOp>::iterator itr = al.CigarData.begin();
        itr != al.CigarData.end(); itr++)
    {
        switch(itr->Type)
        {
            case 'M': // alignment match (can be a sequence match or mismatch)
                iStart += itr->Length;
                break;
            case 'I': // insertion to the reference
                iStart += itr->Length;
                break;
            case 'D': // deletion from the reference
            case 'N':  // skipped region from the reference
                break;
            case 'S':  // soft clipping (clipped sequences present in SEQ)
            case 'H':  // hard clipping (clipped sequences NOT present in SEQ)
                if(enPart == cpAny)
                    return al.QueryBases.substr(iStart, itr->Length);
                else
                {
                    vClispSeq.push_back(al.QueryBases.substr(iStart, itr->Length));
                    iStart += itr->Length;
                }
                break;
            case 'P': // padding (silent deletion from padded reference)
            case '=': // sequence match
            case 'X': // sequence mismatch
                iStart += itr->Length;
                break;
        }
    }
    if(vClispSeq.empty())
        return "";
    else // 在这里我们从来不考虑最长的，我们只考虑该soft clip所在的位置
    {
        if(enPart == cpLeft)
            return *vClispSeq.begin();
        else if(enPart == cpRight)
            return *(vClispSeq.end() - 1);
        else
            return "";
    }
}*/

void WriteFaFile(ofstream& ofs, string& strFilePath, vector<string>& vSeq)
{
    ofs.open(strFilePath.c_str());
    int iIndex = 0;
    for(vector<string>::iterator itr = vSeq.begin(); itr != vSeq.end(); itr++)
    {
        ofs << ">" << IntToStr(iIndex) << endl;
        ofs << *itr << endl;
        iIndex++;
    }
    ofs.close();
}

unsigned long GetFileSize(const char *path)
{
    unsigned long filesize = 0;
    struct stat statbuff;
    if(stat(path, &statbuff) < 0)
    {
        return filesize;
    }
    else
    {
        filesize = statbuff.st_size;
    }
    return filesize;
}
