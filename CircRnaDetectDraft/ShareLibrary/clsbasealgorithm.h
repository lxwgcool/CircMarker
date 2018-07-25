/******************This is head file (*.h) of basical function**************
 * *************************************************************************/

#ifndef CLSBASEALGORITHM_H
#define CLSBASEALGORITHM_H

//#define USEBWA

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;

#ifdef USEBWA
    #include "bamtools/include/api/BamReader.h"
    using namespace BamTools;
#endif


string& ltrim(string &str);
string& rtrim(string &str);
string& trim(string &str);

string IntToStr(int iValue);
string FloatToStr(double fValue, int iDigit=2);

string& ToUpper(string& strValue);
string& ToLower(string& strValue);

string GetReverse(string strOrg);

string GetHigherFolderPath(string strCurPath, int iLevel=1);
string GetCurExeFolderPath();

bool IsMissing(char nt);
char GetComplement(char bp);
string GetReverseCompelement(string strOrg, bool bRevsCompelement = true);

void DisplayString(ofstream& ofs, string& strValue, int iLenPerLine=100);

string::size_type GetNextNPos(string& strOrg, int iCurPos=0, bool bNormalOrientation = true);
bool CheckNPos(char cCurChar);

string GetStringItemValue(string& strValue, string::size_type sztpPos, bool bForward, char cSplit = ' ');

string GetRatio(float fValue, int iPrecision=2, bool bPrintSymbol=true);

void GetAllFilesInFolder(vector<string>& vFileName, const char* czFolderPath);

string GetCompleCode(int iOrg, int iStdDigit=5); // for example: 5 --> 00005; 22 --> 00022

string GetHMSTimeFormat(double dSumSeconds);

bool CheckContainState(int iStart, int iEnd, int iPos,
                       bool bTolerant = false, int iOffSet = 10);

enum En_ClipPart{cpLeft=0, cpRight, cpAny, cpMax};
//string GetSoftClipPart(BamAlignment& al, En_ClipPart enPart=cpAny); // The part that been cut off
//int GetMatchPartLen(BamAlignment& al);

void WriteFaFile(ofstream& ofs, string& strFilePath, vector<string>& vSeq);

unsigned long GetFileSize(const char *path);

#endif // CLSBASEALGORITHM_H
