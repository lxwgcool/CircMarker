#include "clsbwa.h"
#include "clsbasealgorithm.h"
#include "unistd.h" // this is for: get_current_dir_name
#include "stdlib.h"

string ClsBWA::CreateBamBySingleReads(string& strRefPath, string& strReadsPath,
                                      string strBamFileName, string strBamFolder,
                                      bool bMapAll, bool bLooseMatch,
                                      bool bSortByName, int iThreadNum, bool bUseSystemLib)
{
    string strRootPath = ::GetHigherFolderPath(get_current_dir_name());
    string strBWAPath = "bwa";
    if(!bUseSystemLib)
    {
        strBWAPath = "../../ShareLibrary/ThirdPartyTools/bwa/bwa";//strRootPath + "ThirdPartyTools/bwa/bwa";
    }

    if(strBamFolder == "")
        strRootPath += "TempFile/";
    else
    {
        strRootPath += strBamFolder;
        if(strBamFileName != "")  // 如果之前已经做过了，那么我们据不再做了
        {
            string stBamfileName = strRootPath + strBamFileName + ".sorted.bam";
            //Check if the bam file "WholeScaffoldAlign.sorted.bam" has been existed
            if(::access(stBamfileName.c_str(), 0) == 0) // Means it is exsited
            {
                return stBamfileName;
            }
        }
    }
    //2: Build Index for this Reference Fa File
    string strCmd = "";
    string strSpcialFileOfRefIndex = strRefPath + ".bwt";
    if(::access(strSpcialFileOfRefIndex.c_str(), 0) != 0) //This means such file DO NOT exsited
    {
        strCmd = strBWAPath + " index -a bwtsw " + strRefPath;
        system(strCmd.c_str());
    }

    //4: Generate SAM file
    string strSamPath = strRootPath + (strBamFileName == "" ? "Read.sam" : (strBamFileName + ".sam"));
    string strThread = "-t " + IntToStr(iThreadNum) + " ";
    if(bMapAll) //Use map all
    {
        if(bLooseMatch)
            strCmd = strBWAPath + " mem -a -D 0 " + strThread +
                     strRefPath + " " + strReadsPath + " > " + strSamPath;
        else
            strCmd = strBWAPath + " mem -a " + strThread +
                     strRefPath + " " + strReadsPath + " > " + strSamPath;
    }
    else
    {
        if(bLooseMatch)
            strCmd = strBWAPath + " mem -D 0 " + strThread +
                     strRefPath + " " + strReadsPath + " > " + strSamPath;
        else
            strCmd = strBWAPath + " mem " + strThread +
                     strRefPath + " " + strReadsPath + " > " + strSamPath;
    }
    system(strCmd.c_str());

    //5: transfer sam to bam
    string strBamPath = strRootPath + (strBamFileName == "" ? "Read.bam" : (strBamFileName + ".bam"));
    strCmd = "samtools view -bS " + strSamPath + " > " + strBamPath;
    system(strCmd.c_str());

    //6: Sort bam file
    string strSortedBamPath = strRootPath +
                              (strBamFileName == "" ? "Read.sorted.bam" :
                                                      (strBamFileName + ".sorted.bam"));
    //Normal sorted
    string strBamTools = "bamtools";
    if(!bUseSystemLib)
    {
        strBamTools = "../../ShareLibrary/bamtools/bin/bamtools";
    }

    if(bSortByName)
        strCmd = strBamTools + " sort -in " + strBamPath + " -out " +
                 strSortedBamPath + " -byname";
    else
        strCmd = strBamTools + " sort -in " + strBamPath + " -out " + strSortedBamPath;
    system(strCmd.c_str());

    //7:build index file for bam file
    strCmd = strBamTools + " index -in " + strSortedBamPath;
    system(strCmd.c_str());

    //8: delete the bamfile and samfile since they are too large to be stored
    //a) Delete Sam File
    strCmd = "rm -f " + strSamPath;
    system(strCmd.c_str());
    //b) Delete Bam File
    strCmd = "rm -f " + strBamPath;
    system(strCmd.c_str());

    return strSortedBamPath;
}

//do not need to add suffix of Bam_File_Name_Customize
string ClsBWA::CreateBamByPEReads(string& strRefPath, string& strReads1Path, string& strReads2Path,
                                  bool bSortByName, string strBamFolder,
                                  bool bMapAll, bool bLooseMatch,
                                  string strBamFileNameCustmize, int iThreadNum, bool bUseSystemLib)
{
    string strCmd = "";
    string strRootPath = GetHigherFolderPath(get_current_dir_name());

    cout << "-------->" << endl;
    cout << GetHigherFolderPath(get_current_dir_name()) << endl;
    cout << "<--------" << endl;

    string strBWAPath = "bwa";
    if(!bUseSystemLib)
    {
        strBWAPath = "../../ShareLibrary/ThirdPartyTools/bwa/bwa";//strRootPath + "ThirdPartyTools/bwa/bwa";
    }

    if(strBamFolder != "")
    {
        strRootPath += strBamFolder;
        if(strBamFileNameCustmize != "")
        {
            string stBamfileName = strRootPath + strBamFileNameCustmize + ".sorted.bam";
            //Check if the bam file "WholeScaffoldAlign.sorted.bam" has been existed
            if(::access(stBamfileName.c_str(), 0) == 0) // Means it is exsited
            {
                return stBamfileName;
            }
        }
    }
    else
        strRootPath += "TempFile/";

    //Build Index for this Reference Fa File
    string strSpcialFileOfRefIndex = strRefPath + ".bwt";
    if(::access(strSpcialFileOfRefIndex.c_str(), 0) != 0) //This means such file DO NOT exsited
    {
        strCmd = strBWAPath + " index -a bwtsw " + strRefPath;
        system(strCmd.c_str());
    }

    string strSamPath = strRootPath + "Read.sam";
    if(strBamFileNameCustmize != "")
    {
        strSamPath = strRootPath + strBamFileNameCustmize + ".sam";
    }

    string strThread = "-t " + IntToStr(iThreadNum) + " ";
    if(bMapAll) //Use map all
    {
        if(bLooseMatch)
            strCmd = strBWAPath + " mem -a -D 0 " + strThread +
                     strRefPath + " " + strReads1Path + " " + strReads2Path + " > " + strSamPath;
        else
            strCmd = strBWAPath + " mem -a " + strThread +
                     strRefPath + " " + strReads1Path + " " + strReads2Path + " > " + strSamPath;
    }
    else
    {
        if(bLooseMatch)
            strCmd = strBWAPath + " mem -D 0 " + strThread +
                     strRefPath + " " + strReads1Path + " " + strReads2Path + " > " + strSamPath;
        else
            strCmd = strBWAPath + " mem " + strThread +
                     strRefPath + " " + strReads1Path + " " + strReads2Path + " > " + strSamPath;
    }
    cout << endl << "------" << endl;
    cout << strCmd << endl << "------" << endl << endl;

    system(strCmd.c_str());

    //5: transfer sam to bam
    string strBamPath = strRootPath + "Read.bam";
    if(strBamFileNameCustmize != "")
    {
        strBamPath = strRootPath + strBamFileNameCustmize + ".bam";
    }
    strCmd = "samtools view -bS " + strSamPath + " > " + strBamPath;
    system(strCmd.c_str());

    //6: Sort bam file
    string strBamTools = "bamtools";
    string strSortedBamPath = strRootPath + "Read.sorted.bam";
    if(strBamFileNameCustmize != "")
    {
        strSortedBamPath = strRootPath + strBamFileNameCustmize + ".sorted.bam";
    }


    if(!bUseSystemLib)
    {
        strBamTools = "../../ShareLibrary/bamtools/bin/bamtools";
    }

    //check if we need to sort the bam file by name
    if(bSortByName)
        strCmd = strBamTools + " sort -in " +
                strBamPath + " -out " + strSortedBamPath + " -byname";
    else
        strCmd = strBamTools + " sort -in " +
                strBamPath + " -out " + strSortedBamPath;
    system(strCmd.c_str());

    //7:build index file for bam file
    strCmd = strBamTools + " index -in " + strSortedBamPath;
    system(strCmd.c_str());

    //8:Remove the sam file and the bam file
    strCmd = "rm -f " + strSamPath;
    system(strCmd.c_str());
    //b) Delete Bam File
    strCmd = "rm -f " + strBamPath;
    system(strCmd.c_str());

    return strSortedBamPath;
}
